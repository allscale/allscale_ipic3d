#pragma once

#include <iostream>
#include <atomic>
#include <mpi.h>

#include <allscale/utils/serializer.h>
#include <allscale/utils/serializer/vectors.h>
#include <allscale/utils/assert.h>
#include <allscale/utils/vector.h>
#include <allscale/utils/printer/set.h>

#include "ipic3d/app/transfer_buffer.h"

namespace ipic3d {

	namespace partitioning {


		namespace detail {

			std::size_t flattenCellCoordinates(const allscale::utils::Vector<int64_t,3>& size, const int64_t x, const int64_t y, const int64_t z) {
				return (x*size.y + y)*size.z + z;
			}

			// --- Binary splitting distribution ---

			inline int minDivisior(int x) {
				for(int d=2; ; d++) {
					if (x%d == 0) return d;
				}
				return x;
			}

			inline void initDistribution(std::vector<int>& distribution, const allscale::utils::Vector<int64_t,3>& fullSize, int level, int beginRank, int endRank, const coordinate_type& min, const coordinate_type& max) {

				// if only one rank left => base case
				if (beginRank + 1 == endRank) {

					int counter = 0;

					// fill with minRank
					for(int x=min.x; x<max.x; x++) {
						for(int y=min.y; y<max.y; y++) {
							for(int z=min.z; z<max.z; z++) {
								auto pos = flattenCellCoordinates(fullSize,x, y, z);
								distribution[pos] = beginRank;
								counter++;
							}
						}
					}

					// Distribution debugging:
	//				if (isMaster()) std::cout << "Assigning " << min << " - " << max << " to rank " << beginRank << " - total of " << counter << " cells\n";

					return;
				}

				// need to split further
				int splitDim = level % 3;

				// compute number of fragments
				int numRanks = endRank - beginRank;
				int fragments = minDivisior(numRanks);

				int l = max[splitDim] - min[splitDim];
				int share = l / fragments;
				int excess = l % fragments;

				int i = 0;
				int curBegin = min[splitDim];
				int curEnd = curBegin + share + ((i < excess) ? 1 : 0);
				while(curBegin < max[splitDim]) {

					auto curMin = min;
					auto curMax = max;
					curMin[splitDim] = curBegin;
					curMax[splitDim] = curEnd;

					auto curBeginRank = beginRank + i*(numRanks/fragments);
					auto curEndRank = curBeginRank + (numRanks/fragments);
					initDistribution(distribution,fullSize,level+1,curBeginRank, curEndRank, curMin, curMax);

					// compute next sub-fragment of region
					i++;
					curBegin = curEnd;
					curEnd = curBegin + share + ((i < excess) ? 1 : 0);
				}

			}

		} // end namespace detail

		inline std::vector<int> partitionSpaceBinary(int ranks, const allscale::utils::Vector<int64_t,3>& size) {
			auto res = std::vector<int>(size.x * size.y * size.z);

			// initialize recursive
			detail::initDistribution(res,size,0,0,ranks,0,size);

			// done
			return res;
		}

		// --- Z Curve distribution ---

		inline int interleave(int x, int y, int z) {
			int res = 0;
			for(int i=7; i>=0; i--) {
				res = res | ((x >> i) & 0x1);
				res = res << 1;

				res = res | ((y >> i) & 0x1);
				res = res << 1;

				res = res | ((z >> i) & 0x1);
				res = res << 1;
			}
			return res;
		}

		inline std::vector<int> partitionSpaceZCurve(int ranks, const allscale::utils::Vector<int64_t,3>& size) {
			auto res = std::vector<int>(size.x * size.y * size.z);

			using Point = allscale::utils::Vector<int64_t,3>;

			assert_lt(size.x,256);
			assert_lt(size.y,256);
			assert_lt(size.z,256);

			// assign z-curve values
			std::vector<std::pair<int,Point>> points(res.size());
			for(int x = 0; x<size.x; x++) {
				for(int y = 0; y<size.y; y++) {
					for(int z = 0; z<size.z; z++) {
						int val = interleave(x,y,z);
						points[detail::flattenCellCoordinates(size,x,y,z)] = std::make_pair(val,Point(x,y,z));
					}
				}
			}

			// sort points by value
			std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

			// partition
			int fragment = res.size() / ranks;
			int excess = res.size() % ranks;

			int i = 0;
			std::size_t curBegin = 0;
			std::size_t curEnd = curBegin + fragment + (i < excess ? 1 : 0);
			while(curBegin < res.size()) {
				for(std::size_t j=curBegin; j<curEnd; j++) {
					auto pos = points[j].second;
					res[detail::flattenCellCoordinates(size,pos.x,pos.y,pos.z)] = i;
				}

				i++;
				curBegin = curEnd;
				curEnd = curBegin + fragment + (i < excess ? 1 : 0);
			}

			// done
			return res;
		}


		inline std::vector<int> partitionSpace(int ranks, const allscale::utils::Vector<int64_t,3>& size) {
//			return partitionSpaceBinary(ranks,size);
			return partitionSpaceZCurve(ranks,size);
		}
	}

	/**
	 * A global management context of MPI related information
	 * within ipic3d instances.
	 */
	class MPI_Context {

	public:

		using grid_size_t = allscale::utils::Vector<int64_t,3>;

	private:

		static MPI_Context*& getInstancePtr() {
			static MPI_Context* instance = nullptr;
			return instance;
		}

		// the rank of this node
		int rank;

		// the number of ranks
		int size;

		// the distribution of cells
		grid_size_t grid_size;
		std::vector<int> distribution;

		MPI_Context(int argc, char** argv) : rank(0) {
			getInstancePtr() = this;
			MPI_Init(&argc,&argv);
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			MPI_Comm_size(MPI_COMM_WORLD,&size);
			// std::cout << "Initialized MPI on rank " << rank << " / " << size << " ..\n";
		}

	public:

		MPI_Context(const MPI_Context&) = delete;

		MPI_Context(MPI_Context&& other) : rank(other.rank) {
			// move ownership
			getInstancePtr() = this;
		}

		~MPI_Context() {
			if (getInstancePtr() != this) return;
			// std::cout << "Finalizing MPI on rank " << rank << " ..\n";
			MPI_Finalize();
		}

		static MPI_Context init(int argc, char** argv) {
			return MPI_Context(argc,argv);
		}

		static MPI_Context& getInstance() {
			assert_true(getInstancePtr());
			return *getInstancePtr();
		}

		static bool isMaster() {
			return getInstance().rank == 0;
		}

		int getRank() const {
			return rank;
		}

		int getSize() const {
			return size;
		}

		static std::size_t flattenCellCoordinates(const int64_t x, const int64_t y, const int64_t z) {
			auto& size = getInstance().grid_size;
			return (x*size.y + y)*size.z + z;
		}

		static std::size_t flattenCellCoordinates(const grid_size_t& size) {
			return flattenCellCoordinates(size.x, size.y, size.z);
		}

		static void initCellDistribution(const grid_size_t& size) {

			auto& instance = getInstance();
			auto& grid_size = instance.grid_size;
			auto& distribution = instance.distribution;
			auto& group_size = instance.size;

			// check that there has not been an initalization before
			assert_true(distribution.empty()) << "Must only be initialized once!\n";

			// save size
			grid_size = size;

			// create the distribution 'cube'
			distribution = partitioning::partitionSpace(group_size,size);

		}

		static int getRankOf(const grid_size_t& pos) {
			auto& instance = getInstance();
			auto& distribution = instance.distribution;
			return distribution[flattenCellCoordinates(pos)];
		}

		template<typename Body, bool parallel = false>
		static void forEachLocalCell(const Body& body) {
			auto& instance = getInstance();
			auto& size = instance.grid_size;
			int rank = instance.rank;
			// iterate through cells ..
//			#pragma omp parallel for if(parallel) collapse(3)
			for(int x=0; x<size.x; x++) {
				for(int y=0; y<size.y; y++) {
					for(int z=0; z<size.z; z++) {
						// .. and if it is local ..
						if (getRankOf({x,y,z}) == rank) {
							// .. invoke the operation
							body(grid_size_t{x,y,z});
						}
					}
				}
			}
		}

		template<typename Body>
		static void pforEachLocalCell(const Body& body) {
			forEachLocalCell<Body,true>(body);
		}

		template<typename Body, bool parallel = false>
		static void forEachLocalFieldEntry(const Body& body, int additionalEntriesPerDimension) {
			auto& instance = getInstance();
			auto& size = instance.grid_size;
			// iterate through cells ..
//			#pragma omp parallel for if(parallel) collapse(3)
			for(int i=0; i<size.x+additionalEntriesPerDimension; i++) {
				for(int j=0; j<size.y+additionalEntriesPerDimension; j++) {
					for(int k=0; k<size.z+additionalEntriesPerDimension; k++) {
						// .. invoke the operation (all field elements are currently local)
						body(grid_size_t{i,j,k});
					}
				}
			}
		}

		template<typename Body>
		static void pforEachLocalFieldEntry(const Body& body) {
			forEachLocalCell<Body,true>(body);
		}

		static void exchangeBuffers(TransferBuffers& buffers) {

			auto& instance = getInstance();
			auto& size = instance.grid_size;
			auto rank = instance.getRank();

			// sort all particles to export into the list of neighbors

			// mapping from rank to target cells containing the particles to transfer
			std::map<int, std::vector<std::pair<coordinate_type, const std::vector<Particle>*>>> transfers;

			forEachLocalCell([&](auto pos) {
				TransferDirection::forEach([&](auto direction){
					// index particles to be exported
					auto neighborPosition = direction.step(pos, size);
					int targetRank = getRankOf(neighborPosition);
					if(targetRank != rank) {
						transfers[targetRank];	// creates at least an empty list for this target rank
						const auto& particles = buffers.getBuffer(neighborPosition, direction.inverse());
						if (!particles.empty()) {
							transfers[targetRank].push_back({neighborPosition,&particles});
						}
						// clear old buffer state
						auto& buffer = buffers.getBuffer(pos,direction);
						if (!buffer.empty()) buffer.clear();
					}
				});
			});

			std::map<int, allscale::utils::ArchiveWriter> serializers;
			for(const auto& transfer : transfers) {
				auto& out = serializers[transfer.first];
				out.write<std::size_t>(transfer.second.size());
				for(const auto& entry : transfer.second) {
					assert_false(entry.second->empty());
					out.write<coordinate_type>(entry.first);
					out.write<std::vector<Particle>>(*entry.second);
				}
			}

			// send data to all our neighbors - nonblocking
			MPI_Request sendRequests[transfers.size()];
			std::vector<allscale::utils::Archive> archives;
			int sendIndex = 0;
			for(auto& serializer : serializers) {
				archives.emplace_back(std::move(serializer.second).toArchive());
				auto& archive = archives.back();
				auto& buffer = archive.getBuffer();
				MPI_Isend(const_cast<char*>(&buffer[0]), buffer.size(), MPI_BYTE, serializer.first, 0, MPI_COMM_WORLD, &sendRequests[sendIndex++]);
			}

			// receive a message from all our neighbors
			for(std::size_t i = 0; i < transfers.size(); ++i) {
				MPI_Status status;

				// M-version of API not supported by intel MPI
				int flag = 0;
				while (!flag) {
					MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
				}

				int size;
				MPI_Get_count(&status, MPI_BYTE, &size);
				std::vector<char> buffer(size);
				MPI_Recv(&buffer[0], size, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				allscale::utils::Archive archive(buffer);
				allscale::utils::ArchiveReader reader(archive);

				std::size_t numberOfCells = reader.read<std::size_t>();
				for(std::size_t j = 0; j < numberOfCells; ++j) {
					auto cellPosition = reader.read<coordinate_type>();
					auto particles = reader.read<std::vector<Particle>>();
					assert_false(particles.empty());
					assert_eq(instance.getRankOf(cellPosition), rank) << "For cell position: " << cellPosition;
					auto& targetVector = buffers.getBuffer(cellPosition, TransferDirection(0, 0, 0));
					targetVector.insert(targetVector.end(), particles.cbegin(), particles.cend());
				}
			}

			// wait for all sent messages to complete
			MPI_Waitall(transfers.size(), sendRequests, MPI_STATUSES_IGNORE);
		}

	};

} // end namespace ipic3d
