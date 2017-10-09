#pragma once

#include "allscale/api/core/io.h"

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

	// write the output header
	template<typename StreamObject>
	void writeOutputHeader(StreamObject& streamObject) {
		streamObject << "Cycle \t Total Moment \t E energy \t B energy \t Total KE \n";
	}
	
	// write output depending on the set frequency
	template<typename StreamObject>
	void writeOutputData(const int cycle, const int numSteps, Universe& universe, StreamObject& streamObject) {
		if ( (cycle % universe.properties.FieldOutputCycle == 0) || (cycle+1 == numSteps) ) {
			auto getE = [](const auto& field, const auto& index) { return field[index].E; };
			auto getB = [](const auto& field, const auto& index) { return field[index].B; };

			double Eenergy = getEnergy(universe.field, universe.properties, getE);
			double Benergy = getEnergy(universe.field, universe.properties, getB);
			
			// compute total particles kinetic energy
			struct TotalParticleEnergies {
				double total_ke = 0.0;
				double total_mom = 0.0;
			};

			auto map = [&](const coordinate_type& index, TotalParticleEnergies& sum) {
				sum.total_ke += getParticlesKineticEnergy(universe.cells[index]);
				sum.total_mom += getParticlesMomentum(universe.cells[index]);
			};

			auto reduce = [&](const TotalParticleEnergies& a, const TotalParticleEnergies& b) { 
				return TotalParticleEnergies{a.total_ke + b.total_ke, a.total_mom + b.total_mom};
			};
			auto init = []() { return TotalParticleEnergies{0.0,0.0}; };

			auto totalParticleEnergies = allscale::api::user::preduce(coordinate_type(0), universe.cells.size(), map, reduce, init).get();

			streamObject 
				<< cycle << "\t" 
				<< totalParticleEnergies.total_mom << "\t" 
				<< Eenergy << "\t" 
				<< Benergy << "\t" 
				<< totalParticleEnergies.total_ke 
				<< "\n";
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
		
		// create the output file
		auto& manager = allscale::api::core::FileIOManager::getInstance();
		// define the output file name
		std::string outputFilename = universe.properties.outputFileBaseName + "ConservedQuantities.out";
		// create the result file
		auto logFile = manager.createEntry(outputFilename);
		auto outtxt = manager.openOutputStream(logFile);
		writeOutputHeader(outtxt);

		// -- run the simulation --

		// run time loop for the simulation
		for(unsigned i = 0; i < numSteps; ++i) {

			using namespace allscale::api::user;

			// write output to a file: total energy, momentum, E and B total energy
			writeOutputData(i, numSteps, universe, outtxt);

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

		// close the IO manager 
		manager.close(outtxt);

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
