#include <cstdlib>
#include <iostream>

#include "allscale/api/user/data/vector.h"

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
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <config-file>" << std::endl;
		return EXIT_FAILURE;
	}

	// load input configuration
	std::cout << "Loading configuration file \"" << argv[1] << "\" ..." << std::endl;
	auto params = Parameters::read(argv[1]);

	params.dx = 10;
	params.dy = 10;
	params.dz = 10;

	// ----- initialize simulation environment ------

	// setup simulation
	std::cout << "Initializing simulation state ..." << std::endl;

	// initialize universe
	auto universe = createUniverseFromParams(params);

	// get the number of particles in all cells before the simulation begins
	// TODO: disable it for performance runs
	int start_particles = countParticlesInDomain(universe);

	// ----- run the simulation ------

	std::cout << "Running simulation..." << std::endl;

	// -- run the simulation --

	simulateSteps(params.ncycles, universe);

	// ----- finish ------

	// get the number of particles in all cells at the end of the simulation
	// TODO: disable it for performance runs
	int end_particles = countParticlesInDomain(universe);
	if (start_particles != end_particles) {
		std::cout << "[Error]: Periodic boundary conditions on particles were not preserved!\n";
	}

	// be done
	std::cout << "Simulation finished successfully!" << std::endl;
	return EXIT_SUCCESS;
}

