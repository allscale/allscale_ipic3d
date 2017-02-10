#include "ipic3d/simulator.h"

namespace ipic3d {

	using namespace allscale::api::user::data;

	void solveMaxwellEquations(Field&, const Density&);

	void simulateStep(double dt, Cells& cells, Field& field) {
		simulateSteps(1,dt,cells,field);
	}

	void simulateSteps(int numSteps, double dt, Cells& cells, Field& field) {

		// extract size of grid
		auto zero = utils::Coordinate<3>(0);
		auto size = cells.size();
		auto fieldSize = field.size();


		// -- auxiliary structures for communication --

		// the 3-D density field
		Density density(field.size());

		// a 3-D structure collecting contributions of cells to the density
		Grid<DensityCell,3> densityContributions(fieldSize*2);

		// create a buffer for particle transfers
		Grid<std::vector<Particle>,3> particleTransfers(size * 3);	// a grid of buffers for transferring particles between cells


		// -- run the simulation --

		// run time loop for the simulation
		for(int i = 0; i<numSteps; ++i) {

			using namespace allscale::api::user;

			// STEP 1: collect particle contributions

			// project particles to density field
			pfor(zero,size,[&](const utils::Coordinate<3>& pos){
				cells[pos].projectToDensityField(pos,densityContributions);
			});

			// update density field
			pfor(zero,fieldSize,[&](const utils::Coordinate<3>& pos) {
				DensityCell& entry = density[pos];

				// reset density field
				entry.rho = 0.0;
				entry.J = 0.0;

				// aggregate contributions of adjacent cell
				for(int i=0; i<2; i++) {
					for(int j=0; j<2; j++) {
						for(int k=0; k<2; k++) {
							auto& cur = densityContributions[(pos * 2) + utils::Coordinate<3>{i,j,k}];
							entry.rho += cur.rho;
							entry.J += cur.J;
						}
					}
				}
			});


			// STEP 2: solve field equations
			solveMaxwellEquations(field, density);

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 3: project forces to particles and move particles
			pfor(zero,size,[&](const utils::Coordinate<3>& pos){
				cells[pos].moveParticles(pos,field,particleTransfers, dt);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 4: import particles into destination cells
			pfor(zero,size,[&](const utils::Coordinate<3>& pos){
				cells[pos].importParticles(pos,particleTransfers);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

		}

	}


	void solveMaxwellEquations(Field&, const Density&) {
		// TODO: implement
	}


} // end namespace ipic3d
