#pragma once

#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/universe.h"

namespace ipic3d {



	// -------------------------------------------------------------------------------------
	//										Declarations
	// -------------------------------------------------------------------------------------


	template<typename T>
	using Grid = allscale::api::user::data::Grid<T,3>;

	namespace detail {

		struct default_particle_to_field_projector;

		struct default_field_solver;

		struct boris_mover;
	}


	template<
		typename ParticleToFieldProjector 	= detail::default_particle_to_field_projector,
		typename FieldSolver 				= detail::default_field_solver,
		typename ParticleMover 				= detail::boris_mover
	>
	void simulateSteps(const UseCase& useCase, int numSteps, Universe& universe);


	template<
		typename ParticleToFieldProjector 	= detail::default_particle_to_field_projector,
		typename FieldSolver 				= detail::default_field_solver,
		typename ParticleMover 				= detail::boris_mover
	>
	void simulateStep(const UseCase& useCase, Universe& universe) {
		simulateSteps<ParticleToFieldProjector,FieldSolver,ParticleMover>(useCase, 1, universe);
	}




	// -------------------------------------------------------------------------------------
	//										Definitions
	// -------------------------------------------------------------------------------------

	template<
		typename ParticleToFieldProjector,
		typename FieldSolver,
		typename ParticleMover
	>
	void simulateSteps(const UseCase& useCase, int numSteps, Universe& universe) {

		// instantiate operators
		auto particletoFieldProjector = ParticleToFieldProjector();
		auto fieldSolver = FieldSolver();
		auto particleMover = ParticleMover();

		// -- setup simulation --

		// extract size of grid
		auto zero = utils::Coordinate<3>(0);
		auto size = universe.cells.size();
		auto fieldSize = universe.field.size();


		// -- auxiliary structures for communication --

		// the 3-D density field
		Density density(universe.field.size());

		// a 3-D structure collecting contributions of cells to the density
		Grid<DensityCell> densityContributions(fieldSize*2);

		// create a buffer for particle transfers
		Grid<std::vector<Particle>> particleTransfers(size * 3);	// a grid of buffers for transferring particles between cells


		// -- run the simulation --

		// run time loop for the simulation
		for(int i = 0; i<numSteps; ++i) {

			using namespace allscale::api::user;

			// STEP 1: collect particle contributions

			// project particles to density field
			pfor(zero,size,[&](const utils::Coordinate<3>& pos){
				particletoFieldProjector(universe.properties,universe.cells[pos],pos,densityContributions);
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
			// TODO: fieldSolver(universe.field,density,universe.cells);
			pfor(zero,size,[&](const utils::Coordinate<3>& pos){
				fieldSolver(useCase, universe.properties, pos, universe.field);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 3: project forces to particles and move particles
			pfor(zero,size,[&](const utils::Coordinate<3>& pos){
				particleMover(universe.properties,universe.cells[pos],pos,universe.field,particleTransfers);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 4: import particles into destination cells
			pfor(zero,size,[&](const utils::Coordinate<3>& pos){
				importParticles(universe.properties,universe.cells[pos],pos,particleTransfers);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

		}

	}

	namespace detail {

		struct default_particle_to_field_projector {
			void operator()(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, Grid<DensityCell>& densityContributions) const {
				// TODO: make this a free function:
				cell.projectToDensityField(universeProperties,pos,densityContributions);
			}
		};

		struct default_field_solver {
			void operator()(const UseCase& useCase, const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, Field& field) const {
				FieldSolver(useCase, universeProperties, pos, field);
			}

		};

		struct default_particle_mover {
			void operator()(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field, Grid<std::vector<Particle>>& particleTransfers) const {
				moveParticlesFirstOrder(universeProperties,cell,pos,field);
				exportParticles(universeProperties,cell,pos,particleTransfers);
			}
		};

		struct boris_mover {
			void operator()(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field, Grid<std::vector<Particle>>& particleTransfers) const {
				moveParticlesBorisStyle(universeProperties,cell,pos,field);
				exportParticles(universeProperties,cell,pos,particleTransfers);
			}
		};
	}



} // end namespace ipic3d
