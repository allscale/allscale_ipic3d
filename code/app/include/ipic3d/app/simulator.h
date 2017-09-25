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

		//struct forward_field_solver;
		struct leapfrog_field_solver;

		struct boris_mover;
	}


	template<
		typename ParticleToFieldProjector 	= detail::default_particle_to_field_projector,
		typename FieldSolver 				= detail::leapfrog_field_solver,
		typename ParticleMover 				= detail::boris_mover
	>
	void simulateSteps(unsigned numSteps, Universe& universe);


	template<
		typename ParticleToFieldProjector 	= detail::default_particle_to_field_projector,
		typename FieldSolver 				= detail::leapfrog_field_solver,
		typename ParticleMover 				= detail::boris_mover
	>
	void simulateStep(Universe& universe) {
		simulateSteps<ParticleToFieldProjector,FieldSolver,ParticleMover>(1, universe);
	}




	// -------------------------------------------------------------------------------------
	//										Definitions
	// -------------------------------------------------------------------------------------

	// write output depending on the set frequency
	void writeOutput(const int cycle, const int numSteps, Universe& universe) {
		if ( (cycle % universe.properties.FieldOutputCycle == 0) || (cycle+1 == numSteps) ) {
			std::cout << cycle << ' ';
	
			double Eenergy = getEenergy(universe.field, universe.properties);
			double Benergy = getBenergy(universe.field, universe.properties);
	
			// compute total particles kinetic energy
			double total_ke = 0.0;
			double total_mom = 0.0;
			auto zero = utils::Coordinate<3>(0);
			auto size = universe.cells.size();
			allscale::api::user::pfor(zero, size, [&](const utils::Coordinate<3>& pos) {
				total_ke += getParticlesKineticEnergy(universe.cells[pos]);	
				total_mom += getParticlesMomentum(universe.cells[pos]);	
			});

			std::cout << total_mom << ' ';
			std::cout << Eenergy << ' ';
			std::cout << Benergy << ' ';
			std::cout << total_ke << '\n';
		}
	}

	template<
		typename ParticleToFieldProjector,
		typename FieldSolver,
		typename ParticleMover
	>
	void simulateSteps(unsigned numSteps, Universe& universe) {

		// instantiate operators
		auto particleToFieldProjector = ParticleToFieldProjector();
		auto fieldSolver = FieldSolver();
		auto particleMover = ParticleMover();

		// -- setup simulation --

		// extract size of grid
		auto zero = utils::Coordinate<3>(0);
		auto size = universe.cells.size();
		auto densitySize = size + utils::Coordinate<3>(1); // J is defined on nodes
		auto fieldSize = universe.field.size();
		auto fieldStart = utils::Coordinate<3>(1);
		auto fieldEnd = fieldSize - utils::Coordinate<3>(1); // one because a shift due to the boundary conditions


		// -- auxiliary structures for communication --

		// the 3-D density field
		DensityNodes density(densitySize);

		// create a buffer for particle transfers
		Grid<std::vector<Particle>> particleTransfers(size * 3);	// a grid of buffers for transferring particles between cells


		// -- run the simulation --

		// run time loop for the simulation
		for(unsigned i = 0; i < numSteps; ++i) {

			using namespace allscale::api::user;

			// write output to a file: total energy, momentum, E and B total energy
			writeOutput(i, numSteps, universe);

			// STEP 1: collect particle contributions
			// project particles to density field
			pfor(zero, size, [&](const utils::Coordinate<3>& pos) {
				// TODO: this can be improved by adding rho
				// 	J is defined on nodes
				particleToFieldProjector(universe.properties, universe.cells[pos], pos, density);
			});

			// STEP 2: solve field equations
			// update boundaries
			updateFieldsOnBoundaries(universe.field, universe.bcfield);
			
			// TODO: can we call it like that fieldSolver(universe.field,density,universe.cells);
			pfor(fieldStart, fieldEnd, [&](const utils::Coordinate<3>& pos){
				fieldSolver(universe.properties, pos, density, universe.field, universe.bcfield);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 3: project forces to particles and move particles
			pfor(zero, size, [&](const utils::Coordinate<3>& pos){
				particleMover(universe.properties, universe.cells[pos], pos, universe.field, particleTransfers);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 4: import particles into destination cells
			pfor(zero,size,[&](const utils::Coordinate<3>& pos){
				importParticles(universe.properties, universe.cells[pos], pos, particleTransfers);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

		}

	}

	namespace detail {

		struct default_particle_to_field_projector {
			void operator()(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, DensityNodes& density) const {
				projectToDensityField(universeProperties,cell,pos,density);
			}
		};

		struct default_field_solver {
			void operator()(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, DensityNodes& /*density*/, Field& field, BcField& /*bcfield*/) const {
				solveFieldStatically(universeProperties, pos, field);
			}
		};

		struct forward_field_solver {
			void operator()(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, DensityNodes& density, Field& field, BcField& bcfield) const {
				solveFieldForward(universeProperties, pos, density, field, bcfield);
			}
		};

		struct leapfrog_field_solver {
			void operator()(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, DensityNodes& density, Field& field, BcField& bcfield) const {
				solveFieldLeapfrog(universeProperties, pos, density, field, bcfield);
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
