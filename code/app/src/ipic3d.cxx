#include <cstdlib>
#include <iostream>
#include <chrono>

#include "allscale/utils/vector.h"

#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/parameters.h"
#include "ipic3d/app/simulator.h"
#include "ipic3d/app/universe.h"
#include "ipic3d/app/common.h"

#include "ipic3d/app/utils/points.h"

using namespace ipic3d;


namespace {

	using namespace allscale::utils;

	const uint64_t NUM_TIME_STEPS = 100;
	const double DELTA_T = 0.15;

	const Vector3<int64_t> GRID_SIZE { 16, 16, 16 };

	const Vector3<double> CELL_WIDTH { 10, 10, 10 };

	const Vector3<double> UNIVERSE_SIZE {
		GRID_SIZE.x * CELL_WIDTH.x,
		GRID_SIZE.y * CELL_WIDTH.y,
		GRID_SIZE.z * CELL_WIDTH.z
	};

	int processUniverse(Universe& universe, std::uint64_t numParticles, std::uint64_t numTimeSteps) {
		const bool dumpParticles = std::getenv("DUMP_PARTICLE_POSITION");

		// --- sample output initial state --
		if (dumpParticles) outputParticlePositionSamples(universe.cells, std::string("t_begin.txt"), 1000);


		// ----- run the simulation ------
		std::cout << "Running simulation..." << std::endl;

		auto start = std::chrono::high_resolution_clock::now();
		simulateSteps(numTimeSteps, universe);
		auto end = std::chrono::high_resolution_clock::now();

		std::cout << "Simulation Finished" << std::endl;

		double s = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() / 1000.0;
		std::cout << "Simulation took " << s << "s\n";
		std::cout << "Throughput: " << (numTimeSteps * numParticles) / s << " particles/s \n";

		if (dumpParticles) outputParticlePositionSamples(universe.cells, std::string("t_end.txt"), 1000);

		return EXIT_SUCCESS;
	}

	template<typename Distribution>
	int processDistribution(const Distribution& dist, std::uint64_t numParticles, std::uint64_t numTimeSteps) {

		// set up universe
		UniverseProperties universeProperties;
		universeProperties.cellWidth = CELL_WIDTH;
		universeProperties.size = GRID_SIZE;
		universeProperties.dt = DELTA_T;

		// set up init properties
		InitProperties initProperties;
		initProperties.driftVelocity.push_back(0);

		// initialize universe
		auto universe = createUniverseFromDistribution(universeProperties,initProperties,numParticles,dist);

		// -- run the simulation --

		return processUniverse(universe,numParticles,numTimeSteps);

	}

	int processUniform(std::uint64_t numParticles, std::uint64_t numTimeSteps) {
		return processDistribution(
			ipic3d::distribution::uniform<>(
					Vector3<double> { 0 ,0 ,0 },
					UNIVERSE_SIZE,
					Vector3<double> { -0.2 , -0.2 , -0.2},
					Vector3<double> { +0.2 , +0.2 , +0.2}
			),
			numParticles,
			numTimeSteps
		);
	}

	int processCluster(std::uint64_t numParticles, std::uint64_t numTimeSteps) {
		return processDistribution(
			ipic3d::distribution::normal<>(
					UNIVERSE_SIZE/2,
					UNIVERSE_SIZE/5,
					Vector3<double> { -0.2 , -0.2 , -0.2},
					Vector3<double> { +0.2 , +0.2 , +0.2}
			),
			numParticles,
			numTimeSteps
		);
	}

	int processExplosion(std::uint64_t numParticles, std::uint64_t numTimeSteps) {
		return processDistribution(
			ipic3d::distribution::make_spherical(
				ipic3d::distribution::generic_particle_generator<
						ipic3d::distribution::vector::uniform,		// < position distribution
						ipic3d::distribution::vector::normal,		// < velocity distribution
						ipic3d::distribution::species::electron		// < species distribution
					>(
						{ UNIVERSE_SIZE * 0.4, UNIVERSE_SIZE * 0.6, 0 },
						{ Vector3<double> { 0 }, Vector3<double> { 1.5 }, 1 },
						{}
				),
				UNIVERSE_SIZE/2,
				UNIVERSE_SIZE.x/10
			),
			numParticles,
			numTimeSteps
		);
	}

	int processBeam(std::uint64_t numParticles, std::uint64_t numTimeSteps) {
		return processDistribution(
			ipic3d::distribution::normal<>(
					UNIVERSE_SIZE/100,
					UNIVERSE_SIZE/500,
					UNIVERSE_SIZE / numTimeSteps * 0.95 / DELTA_T,
					UNIVERSE_SIZE / numTimeSteps * 0.95 / DELTA_T
			),
			numParticles,
			numTimeSteps
		);
	}

	int processBenchmark(const std::string& config) {

		std::cout << "Processing benchmark " << config << "\n";

		// utility for help messages
		auto showHelp = [](){
			std::cout << "Benchmark flags: -X:N\n"
					  << "      X ... benchmark type:\n"
					  << "           U ... uniform\n"
					  << "           C ... cluster\n"
					  << "           E ... explosion\n"
					  << "           B ... beam\n"
					  << "      N ... total number of particles:\n";
			return EXIT_FAILURE;
		};

		// check the format
		if (config.length() <= 3 || config[0] != '-' || config[2] != ':') {
			return showHelp();
		}

		// all benchmarks comprise 100 time steps
		std::uint64_t numTimeSteps = NUM_TIME_STEPS;

		// parse number of particles
		std::uint64_t numParticles = std::atoll(config.substr(3,config.length()).c_str());

		// dispatch based on use case
		switch(config[1]) {
			case 'U': return processUniform(numParticles,numTimeSteps);
			case 'C': return processCluster(numParticles,numTimeSteps);
			case 'E': return processExplosion(numParticles,numTimeSteps);
			case 'B': return processBeam(numParticles,numTimeSteps);
			default:  return showHelp();
		}
	}
}


int main(int argc, char** argv) {

	// ----- load and parse simulation parameters ------

	// check the passed arguments
	if (argc != 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
		std::cout << "Usage: ./ipic3d <config-file>" << std::endl;
		return EXIT_FAILURE;
	}

	std::string inputFilename = argv[1];

	// add benchmark support
	if (inputFilename[0] == '-') {
		return processBenchmark(inputFilename);
	}


	// load input configuration
	std::cout << "Loading configuration file \"" << inputFilename << "\" ..." << std::endl;
	auto params = Parameters(inputFilename);

	// ----- initialize simulation environment ------

	// setup simulation
	std::cout << "Initializing simulation state ..." << std::endl;

	// remove preceeding path from filename and file suffix, keep only file name itself
	const auto sepPos = inputFilename.find_last_of("/\\");
	std::string baseName = inputFilename.substr(sepPos + 1, inputFilename.find_last_of('.') - sepPos - 1);

	// initialize universe
	auto universe = createUniverseFromParams(params, baseName);

	// dump samples of particle positions for visualization
	outputParticlePositionSamples(universe.cells, "p_pos.txt");

	// get the number of particles in all cells before the simulation begins for error checking
	assert_decl(int start_particles = countParticlesInDomain(universe));

	// ----- run the simulation ------

	std::cout << "Running simulation..." << std::endl;

	// -- run the simulation --

	simulateSteps(params.ncycles, universe);

	// ----- finish ------

	// get the number of particles in all cells at the end of the simulation for error checking
	assert_decl(int end_particles = countParticlesInDomain(universe));
	assert_eq(start_particles, end_particles) << "[Error]: Periodic boundary conditions on particles were not preserved!";

	std::cout << "Simulation finished successfully, producing output data..." << std::endl;

	std::string outputFilename = baseName + ".out";
	auto& manager = allscale::api::core::FileIOManager::getInstance();
	auto text = manager.createEntry(outputFilename);
	auto out = manager.openOutputStream(text);
	outputNumberOfParticlesPerCell(universe.cells, out);
	outputFieldGrids(universe.field, universe.bcfield, out);

	// be done
	return EXIT_SUCCESS;
}

