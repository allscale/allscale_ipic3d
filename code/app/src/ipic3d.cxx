#include <cstdlib>
#include <iostream>

#include "allscale/utils/vector.h"

#include "ipic3d/app/cell.h"
#include "ipic3d/app/parameters.h"
#include "ipic3d/app/simulator.h"
#include "ipic3d/app/universe.h"
#include "ipic3d/app/common.h"

#include "ipic3d/app/utils/points.h"

using namespace ipic3d;

int main(int argc, char** argv) {

	// ----- load and parse simulation parameters ------

	// check the passed arguments
	if (argc != 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
		std::cout << "Usage: ./ipic3d <config-file>" << std::endl;
		return EXIT_FAILURE;
	}

	std::string inputFilename = argv[1];

	// load input configuration
	std::cout << "Loading configuration file \"" << inputFilename << "\" ..." << std::endl;
	auto params = Parameters(inputFilename);

	// ----- initialize simulation environment ------

	// setup simulation
	std::cout << "Initializing simulation state ..." << std::endl;

	// initialize universe
	auto universe = createUniverseFromParams(params);

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

	// remove preceeding path from filename and file suffix, keep only file name itself
	const auto sepPos = inputFilename.find_last_of("/\\");
	std::string baseName = inputFilename.substr(sepPos + 1, inputFilename.find_last_of('.') - sepPos - 1);
	std::string outputFilename = baseName + ".out";
	outputNumberOfParticlesPerCell(universe.cells, outputFilename);

	// be done
	std::cout << "Simulation finished successfully!" << std::endl;
	return EXIT_SUCCESS;
}

