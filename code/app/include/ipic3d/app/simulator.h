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

		struct default_particle_mover;
	}


	template<
		typename ParticleToFieldProjector 	= detail::default_particle_to_field_projector,
		typename FieldSolver 				= detail::default_field_solver,
		typename ParticleMover 				= detail::default_particle_mover
	>
	void simulateSteps(unsigned numSteps, Universe& universe);


	template<
		typename ParticleToFieldProjector 	= detail::default_particle_to_field_projector,
		typename FieldSolver 				= detail::default_field_solver,
		typename ParticleMover 				= detail::default_particle_mover
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
		if ( universe.properties.FieldOutputCycle > 0 && ((cycle % universe.properties.FieldOutputCycle == 0) || (cycle+1 == numSteps)) ) {
			auto getE = [](const auto& field, const auto& index) { return field[index].E; };
			auto getB = [](const auto& field, const auto& index) { return (field[index].B + field[index].Bext); };

			double Eenergy = getFieldEnergy(universe.field, universe.properties, getE);
			double Benergy = getFieldEnergy(universe.field, universe.properties, getB);

			double totalParticlesMomentum = getTotalParticlesEnergy(universe.cells, getParticlesMomentum);
			double totalParticlesKineticEnergy = getTotalParticlesEnergy(universe.cells, getParticlesKineticEnergy);

			streamObject 
				<< cycle << "\t" 
				<< totalParticlesMomentum << "\t" 
				<< Eenergy << "\t" 
				<< Benergy << "\t" 
				<< totalParticlesKineticEnergy 
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
		//auto particleToFieldProjector = ParticleToFieldProjector();
		//auto fieldSolver = FieldSolver();
		auto particleMover = ParticleMover();

		// -- setup simulation --

		// extract size of grid
		auto zero = utils::Coordinate<3>(0);
		auto size = universe.cells.size();
		//auto densitySize = universe.currentDensity.size();
		//auto fieldSize = universe.field.size();
		//auto fieldStart = utils::Coordinate<3>(1);
		//auto fieldEnd = fieldSize - utils::Coordinate<3>(1); // one because a shift due to the boundary conditions


		// -- auxiliary structures for communication --

		// create a buffer for particle transfers
		Grid<std::vector<Particle>> particleTransfers(size * 3);	// a grid of buffers for transferring particles between cells

#ifdef ENABLE_DEBUG_OUTPUT
		// create the output file
		auto& manager = allscale::api::core::FileIOManager::getInstance();
		// define the output file name
		std::string outputFilename = universe.properties.outputFileBaseName + "ConservedQuantities.out";
		// create the result file
		auto logFile = manager.createEntry(outputFilename);
		auto outtxt = manager.openOutputStream(logFile);

		writeOutputHeader(outtxt);
#endif

		// -- run the simulation --

		// run time loop for the simulation
		for(unsigned i = 0; i < numSteps; ++i) {

			using namespace allscale::api::user::algorithm;

#ifdef ENABLE_DEBUG_OUTPUT
			// write output to a file: total energy, momentum, E and B total energy
			writeOutputData(i, numSteps, universe, outtxt);
#endif

			// STEP 1: collect particle contributions
			// project particles to current density
			// use size, not densitySize, because all Js in the cell will be updated
			//pfor(zero, densitySize, [&](const utils::Coordinate<3>& pos) {
			//	// TODO: this can be improved by adding rho
			//	// 	J is defined on nodes
			//	particleToFieldProjector(universe.properties, universe.cells, pos, universe.currentDensity);
			//});

			// STEP 2: solve field equations
			// update boundaries
			updateFieldsOnBoundaries(universe.field, universe.bcfield);

			// TODO: can we call it like that fieldSolver(universe.field,density,universe.cells);
			//pfor(fieldStart, fieldEnd, [&](const utils::Coordinate<3>& pos){
			//	fieldSolver(universe.properties, pos, universe.currentDensity, universe.field, universe.bcfield);
			//});

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 3: project forces to particles and move particles
			pfor(zero, size, [&](const utils::Coordinate<3>& pos){
				particleMover(universe.properties, universe.cells[pos], pos, universe.field, particleTransfers);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 4: import particles into destination cells
			pfor(zero, size, [&](const utils::Coordinate<3>& pos){
				importParticles(universe.properties, universe.cells[pos], pos, particleTransfers);
			});

			// -- implicit global sync - TODO: can this be eliminated? --

		}

#ifdef ENABLE_DEBUG_OUTPUT
		// close the IO manager 
		manager.close(outtxt);
#endif

	}

	namespace detail {

		struct default_particle_to_field_projector {
			void operator()(const UniverseProperties& universeProperties, const Cells& cells, const utils::Coordinate<3>& pos, CurrentDensity& density) const {
				projectToDensityField(universeProperties,cells,pos,density);
			}
		};

		struct default_field_solver {
			void operator()(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, CurrentDensity& /*density*/, Field& field, BcField& /*bcfield*/) const {
				solveFieldStatically(universeProperties, pos, field);
			}
		};

		struct forward_field_solver {
			void operator()(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, CurrentDensity& density, Field& field, BcField& bcfield) const {
				solveFieldForward(universeProperties, pos, density, field, bcfield);
			}
		};

		struct leapfrog_field_solver {
			void operator()(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, CurrentDensity& density, Field& field, BcField& bcfield) const {
				solveFieldLeapfrog(universeProperties, pos, density, field, bcfield);
			}
		};

		struct default_particle_mover {
			void operator()(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field, Grid<std::vector<Particle>>& particleTransfers) const {
				moveParticles(universeProperties,cell,pos,field);
				exportParticles(universeProperties,cell,pos,particleTransfers);
			}
		};
	}

} // end namespace ipic3d
