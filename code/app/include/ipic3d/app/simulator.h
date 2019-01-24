#pragma once

#include <chrono>

#include "allscale/api/core/io.h"

#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/transfer_buffer.h"
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

		struct forward_field_solver;

		struct leapfrog_field_solver;

		struct default_particle_mover;
	}

	struct DurationMeasurement {
		double firstStep;
		double remainingSteps;
	};

	template<
		typename ParticleToFieldProjector 	= detail::default_particle_to_field_projector,
		typename FieldSolver 				= detail::default_field_solver,
		typename ParticleMover 				= detail::default_particle_mover
	>
	DurationMeasurement simulateSteps(std::uint64_t numSteps, Universe& universe);


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
	void writeOutputData(const int cycle, const int numSteps, Universe& universe, StreamObject& streamObject, char* fileName) {
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

		if ( universe.properties.ParticleOutputCycle > 0 && ( (cycle == 0) || ((cycle+1) % universe.properties.ParticleOutputCycle == 0) ) ) {
			int t = (cycle+1) / universe.properties.ParticleOutputCycle;
			std::cout << cycle << ' ' << t << '\n';
			char fileNamet[50];
			std::sprintf(fileNamet, "%s%06d", fileName, t);

			// open file and dump results
			auto out = std::fstream(fileNamet, std::ios_base::out);
			out << "t,x,y,z,density\n";
			for(int x = 0; x < universe.properties.size.x; x++) {
				for(int y = 0; y < universe.properties.size.y; y++) {
					for(int z = 0; z < universe.properties.size.z; z++) {
						double dx = x * universe.properties.cellWidth.x;
						double dy = y * universe.properties.cellWidth.y;
						double dz = z * universe.properties.cellWidth.z;
						
						out << t << "," << dx << "," << dy << "," << dz << "," << universe.cells[coordinate_type{x,y,z}].particles.size() << "\n";
					}
				}
			}
		}
	}

	namespace {

		// temporarily moved to function due to a bug in AllScale compiler
		template <typename T>
		double getTimeCount(T duration) {
			return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() / 1000.0f;
		}

	}

	template<
		typename ParticleToFieldProjector,
		typename FieldSolver,
		typename ParticleMover
	>
	DurationMeasurement simulateSteps(std::uint64_t numSteps, Universe& universe) {

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
		TransferBuffers particleTransfers(size);

		// create a grid of buffers for density projection from particles to grid nodes
		Grid<DensityNode> densityContributions(size * 2);
		
#ifdef ENABLE_DEBUG_OUTPUT
		// create the output file
		auto& manager = allscale::api::core::FileIOManager::getInstance();
		// define the output file name
		std::string outputFilename = universe.properties.outputFileBaseName + "ConservedQuantities.out";
		// create the result file
		auto logFile = manager.createEntry(outputFilename);
		auto outtxt = manager.openOutputStream(logFile);

		writeOutputHeader(outtxt);

		// use the current time to create a unique file name
		auto timeStamp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
		char fileName[50];
		std::sprintf(fileName,"result_%ld.csv.", timeStamp);
#endif

		// -- run the simulation --
		
		auto start = std::chrono::high_resolution_clock::now();
		auto endFirst = start;

		allscale::api::user::algorithm::detail::loop_reference<utils::Coordinate<3>> ref;

		// run time loop for the simulation
		for(std::uint64_t i = 0; i < numSteps; ++i) {

			using namespace allscale::api::user::algorithm;

#ifdef ENABLE_DEBUG_OUTPUT
			// write output to a file: total energy, momentum, E and B total energy
			writeOutputData(i, numSteps, universe, outtxt, fileName);
#endif
			// STEP 1a: collect particle density contributions and store in buffers
			//pfor(zero, size, [&](const utils::Coordinate<3>& pos) {
			//	// TODO: this can be improved by adding rho
			//	// 	J is defined on nodes
			//	particleToFieldProjector(universe.properties, universe.cells[pos], pos, densityContributions);
			//});

			//// STEP 1b: aggregate densities in buffers to density nodes
			//pfor(zero, universe.currentDensity.size(), [&](const utils::Coordinate<3>& pos) {
			//	// 	J is defined on nodes
			//	aggregateDensityContributions(universe.properties, densityContributions, pos, universe.currentDensity[pos]);
			//});

			// STEP 2: solve field equations
			// update boundaries
			//updateFieldsOnBoundaries(universe.field, universe.bcfield);

			// TODO: can we call it like that fieldSolver(universe.field,density,universe.cells);
			//pfor(fieldStart, fieldEnd, [&](const utils::Coordinate<3>& pos){
			//	fieldSolver(universe.properties, pos, universe.currentDensity, universe.field, universe.bcfield);
			//});

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 3: project forces to particles and move particles
			ref = pfor(zero, size, [particleMover,&universe,&particleTransfers](const utils::Coordinate<3>& pos){
				particleMover(universe.properties, universe.cells[pos], pos, universe.field, particleTransfers);
			}, allscale::api::user::algorithm::full_neighborhood_sync(ref));

			// -- implicit global sync - TODO: can this be eliminated? --

			// STEP 4: import particles into destination cells
			ref = pfor(zero, size, [&](const utils::Coordinate<3>& pos){
				importParticles(universe.properties, universe.cells[pos], pos, particleTransfers);
			}, allscale::api::user::algorithm::full_neighborhood_sync(ref));

			// -- implicit global sync - TODO: can this be eliminated? --
			
			if(i == 0) {
				endFirst = std::chrono::high_resolution_clock::now();
			}

		}

		// wait for completion
		ref.wait();

		auto endAll = std::chrono::high_resolution_clock::now();
		auto durationFirst = endFirst - start;
		auto durationRemaining = endAll - endFirst;

#ifdef ENABLE_DEBUG_OUTPUT
		// close the IO manager 
		manager.close(outtxt);
#endif
		return { getTimeCount(durationFirst), getTimeCount(durationRemaining) };
	}

	namespace detail {

		struct default_particle_to_field_projector {
			void operator()(const UniverseProperties& universeProperties, const Cells& cells, const utils::Coordinate<3>& pos, CurrentDensity& densityContributions) const {
				projectToDensityField(universeProperties, cells, pos, densityContributions);
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
			void operator()(const UniverseProperties& properties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field, TransferBuffers& particleTransfers) const {
				moveParticles(properties, cell, pos, field);
				exportParticles(properties, cell, pos, particleTransfers);
			}
		};
	}

} // end namespace ipic3d
