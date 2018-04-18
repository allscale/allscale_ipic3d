#include <cstdlib>
#include <iostream>

#include "ipic3d/app/benchmark.h"
#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/parameters.h"
#include "ipic3d/app/simulator.h"
#include "ipic3d/app/universe.h"

#include "ipic3d/app/mpi_context.h"

using namespace ipic3d;

int main(int argc, char** argv) {

	// startup MPI
	auto context = MPI_Context::init(argc,argv);

	// ----- load and parse simulation parameters ------

	// check the passed arguments
	if (argc != 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
		if (MPI_Context::isMaster()) std::cout << "Usage: ./ipic3d <config-file>" << std::endl;
		return EXIT_SUCCESS;
	}

	std::string inputFilename = argv[1];

	// add benchmark support
	if (inputFilename[0] == ':') {
		processBenchmark(inputFilename);
		return EXIT_SUCCESS;
	}

	// load input configuration
	std::cout << "Loading configuration file \"" << inputFilename << "\" ..." << std::endl;
	auto params = Parameters(inputFilename);

	// ----- initialize simulation environment ------

	// setup simulation
	if (MPI_Context::isMaster()) std::cout << "Initializing simulation state ..." << std::endl;

	// remove preceding path from filename and file suffix, keep only file name itself
	const auto sepPos = inputFilename.find_last_of("/\\");
	std::string baseName = inputFilename.substr(sepPos + 1, inputFilename.find_last_of('.') - sepPos - 1);

	// initialize universe
	auto universe = createUniverseFromParams(params, baseName);

#ifdef ENABLE_DEBUG_OUTPUT
	// get the number of particles in all cells before the simulation begins for error checking
	assert_decl(int start_particles = countParticlesInDomain(universe));
#endif

	if (MPI_Context::isMaster()) std::cout << "Running simulation..." << std::endl;

	// -- run the simulation --

	simulateSteps(params.ncycles, universe);

	// ----- finish ------

#ifdef ENABLE_DEBUG_OUTPUT
	// get the number of particles in all cells at the end of the simulation for error checking
	assert_decl(int end_particles = countParticlesInDomain(universe));
	assert_eq(start_particles, end_particles) << "[Error]: Periodic boundary conditions on particles were not preserved!";
#endif

	if (MPI_Context::isMaster()) std::cout << "Simulation finished successfully, producing output data..." << std::endl;

	std::string outputFilename = baseName + ".out";
	outputNumberOfParticlesPerCell(universe.cells, outputFilename);
	//outputFieldGrids(universe.field, universe.bcfield, outputFilename);

	// be done
	return EXIT_SUCCESS;
}
