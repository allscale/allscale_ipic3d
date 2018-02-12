#include <cstdlib>
#include <iostream>

#include "ipic3d/app/benchmark.h"
#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/parameters.h"
#include "ipic3d/app/simulator.h"
#include "ipic3d/app/universe.h"

using namespace ipic3d;

int main(int argc, char** argv) {

	// ----- load and parse simulation parameters ------

	// check the passed arguments
	if (argc != 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
		std::cout << "Usage: ./ipic3d <config-file>" << std::endl;
		return EXIT_FAILURE;
	}

	std::string inputFilename = argv[1];

	// add benchmark support
	if (inputFilename[0] == ':') {
		return processBenchmark(inputFilename);
	}

	// load input configuration
	std::cout << "Loading configuration file \"" << inputFilename << "\" ..." << std::endl;
	auto params = Parameters(inputFilename);

	// ----- initialize simulation environment ------

	// setup simulation
	std::cout << "Initializing simulation state ..." << std::endl;

	// remove preceding path from filename and file suffix, keep only file name itself
	const auto sepPos = inputFilename.find_last_of("/\\");
	std::string baseName = inputFilename.substr(sepPos + 1, inputFilename.find_last_of('.') - sepPos - 1);

	// initialize universe
	auto universe = createUniverseFromParams(params, baseName);

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
//	outputNumberOfParticlesPerCell(universe.cells, out);
//	outputFieldGrids(universe.field, universe.bcfield, out);

	// be done
	return EXIT_SUCCESS;
}
