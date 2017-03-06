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
	void simulateSteps(unsigned numSteps, Universe& universe);


	template<
		typename ParticleToFieldProjector 	= detail::default_particle_to_field_projector,
		typename FieldSolver 				= detail::default_field_solver,
		typename ParticleMover 				= detail::boris_mover
	>
	void simulateStep(Universe& universe) {
		simulateSteps<ParticleToFieldProjector,FieldSolver,ParticleMover>(1, universe);
	}




	// -------------------------------------------------------------------------------------
	//										Definitions
	// -------------------------------------------------------------------------------------

	template<
		typename ParticleToFieldProjector,
		typename FieldSolver,
		typename ParticleMover
	>
	void simulateSteps(unsigned numSteps, Universe& universe) {

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
		Density density(size);

		// create a buffer for particle transfers
		Grid<std::vector<Particle>> particleTransfers(size * 3);	// a grid of buffers for transferring particles between cells


		// -- run the simulation --

		// run time loop for the simulation
		for(unsigned i = 0; i<numSteps; ++i) {

			using namespace allscale::api::user;

			// STEP 1: collect particle contributions

			// project particles to density field
			pfor(zero,size,[&](const utils::Coordinate<3>& pos) {
				DensityCell& entry = density[pos];

				particletoFieldProjector(universe.properties,universe.cells[pos],pos,entry);
			});


			// STEP 2: solve field equations
			// TODO: fieldSolver(universe.field,density,universe.cells);
			pfor(zero,fieldSize,[&](const utils::Coordinate<3>& pos){
				fieldSolver(universe.properties, pos, universe.field);
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
			void operator()(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& /*pos*/, DensityCell& densityContributions) const {
				projectToDensityField(universeProperties,cell,densityContributions);
			}
		};

		struct default_field_solver {
			void operator()(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, Field& field) const {
				solveFieldStatically(universeProperties, pos, field);
			}

		};

		struct leapfrog_field_solver {
			void operator()(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, Field& field) const {
				solveFieldLeapfrog(universeProperties, pos, field);
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
