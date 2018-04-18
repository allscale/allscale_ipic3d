#pragma once

#include <iostream>
#include <mpi.h>

#include <allscale/utils/assert.h>
#include <allscale/utils/vector.h>

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
			for(int i=0; i<size.x; i++) {
				for(int j=0; j<size.y; j++) {
					for(int k=0; k<size.z; k++) {
						// TODO: find some better distribution
						auto pos = (i*size.y + j)*size.z + k;
						distribution[pos] = pos % group_size;
					}
				}
			}

		}

		template<typename Body, bool parallel = false>
		static void forEachLocalCell(const Body& body) {
			auto& instance = getInstance();
			auto& size = instance.grid_size;
			auto& dist = instance.distribution;
			int rank = instance.rank;
			// iterate through cells ..
//			#pragma omp parallel for if(parallel) collapse(3)
			for(int i=0; i<size.x; i++) {
				for(int j=0; j<size.y; j++) {
					for(int k=0; k<size.z; k++) {
						// .. and if it is local ..
						if (dist[(i*size.y + j)*size.z + k] == rank) {
							// .. invoke the operation
							body(grid_size_t{i,j,k});
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
		static void forEachLocalFieldEntry(const Body& body) {
			auto& instance = getInstance();
			auto& size = instance.grid_size;
			// iterate through cells ..
//			#pragma omp parallel for if(parallel) collapse(3)
			for(int i=0; i<size.x+1; i++) {
				for(int j=0; j<size.y+1; j++) {
					for(int k=0; k<size.z+1; k++) {
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

	};

} // end namespace ipic3d
