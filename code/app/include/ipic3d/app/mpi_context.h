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
			distribution.resize(size.x * size.y * size.z);

			// setup distribution (for now, super stupid)
			for(int x=0; x<size.x; x++) {
				for(int y=0; y<size.y; y++) {
					for(int z=0; z<size.z; z++) {
						// TODO: find some better distribution
						auto pos = flattenCellCoordinates(x, y, z);
						distribution[pos] = pos % group_size;
					}
				}
			}

		}

		static int getRankOf(const grid_size_t& pos) {
			auto& instance = getInstance();
			auto& size = instance.grid_size;
			auto& distribution = instance.distribution;
			return distribution[flattenCellCoordinates(pos)];
		}

		template<typename Body, bool parallel = false>
		static void forEachLocalCell(const Body& body) {
			auto& instance = getInstance();
			auto& size = instance.grid_size;
			auto& dist = instance.distribution;
			int rank = instance.rank;
			// iterate through cells ..
//			#pragma omp parallel for if(parallel) collapse(3)
			for(int x=0; x<size.x; x++) {
				for(int y=0; y<size.y; y++) {
					for(int z=0; z<size.z; z++) {
						// .. and if it is local ..
						if (dist[flattenCellCoordinates(x, y, z)] == rank) {
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
			std::map<int, std::map<coordinate_type, std::vector<Particle>>> transfers;

			forEachLocalCell([&](auto pos) {
				TransferDirection::forEach([&](auto direction){

					// index particles to be exported
					auto neighborPosition = direction.step(pos, size);
					int targetRank = getRankOf(neighborPosition);
					if(targetRank != rank) {
						const auto& particles = buffers.getBuffer(neighborPosition, direction.inverse());
						auto& targetVector = transfers[targetRank][neighborPosition];
						targetVector.insert(targetVector.end(), particles.cbegin(), particles.cend());

						// clear old buffer state
						buffers.getBuffer(pos,direction).clear();
					}


				});
			});

//			if(isMaster()) {
//				for(auto& outer : transfers) {
//					std::cout << "Target Rank: " << outer.first << std::endl;
//					for(auto& inner : outer.second) {
//						std::cout << "  Target Cell: " << inner.first << std::endl;
//						std::cout << "    Particles: " << inner.second << std::endl;
//					}
//				}
//			}

			std::map<int, allscale::utils::ArchiveWriter> serializers;
			for(const auto& transfer : transfers) {
				auto& out = serializers[transfer.first];
				out.write<std::size_t>(transfer.second.size());
				for(const auto& entry : transfer.second) {
					out.write<coordinate_type>(entry.first);
					out.write<std::vector<Particle>>(entry.second);
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
			for(int i = 0; i < transfers.size(); ++i) {
				MPI_Status status;
				MPI_Message message;
				MPI_Mprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &message, &status);
				int size;
				MPI_Get_count(&status, MPI_BYTE, &size);
				std::vector<char> buffer(size);
				MPI_Mrecv(&buffer[0], size, MPI_BYTE, &message, MPI_STATUS_IGNORE);

				allscale::utils::Archive archive(buffer);
				allscale::utils::ArchiveReader reader(archive);

				std::size_t numberOfCells = reader.read<std::size_t>();
				for(std::size_t j = 0; j < numberOfCells; ++j) {
					auto cellPosition = reader.read<coordinate_type>();
					auto particles = reader.read<std::vector<Particle>>();
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
