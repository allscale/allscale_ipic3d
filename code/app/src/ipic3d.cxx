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

	// initialize initial properties
	InitProperties initProperties = InitProperties(params);
	std::cout << initProperties;

	// initialize universe properties
	UniverseProperties universeProperties = UniverseProperties(params);
	universeProperties.outputFileBaseName = baseName;
	universeProperties.dt = 0.01;
	universeProperties.speedOfLight = 299792458;
	//int R = 16;
	//universeProperties.size = { R, R, R };
	universeProperties.planetRadius = 6378137; // meter (Earth radius) 
	universeProperties.cellWidth *= universeProperties.planetRadius;
	//universeProperties.FieldOutputCycle = 0;

	// these parameters are required for computations
	//universeProperties.ParticleOutputCycle = 10;

	auto universeSize = elementwiseProduct(universeProperties.cellWidth, universeProperties.size);
	universeProperties.objectCenter = { 0.0, 0.0, 0.0 };
	universeProperties.origin = universeProperties.objectCenter - universeSize / 2.0; 
	universeProperties.externalMagneticField = { 0.0, 0.0, 3.07e-5 };
	universeProperties.useCase = UseCase::Dipole;

	std::cout << universeProperties;

	// initialize universe
	int numParticles = universeProperties.size.x * universeProperties.size.y * universeProperties.size.z;
	numParticles = numParticles * initProperties.particlesPerCell[0].x * initProperties.particlesPerCell[0].y * initProperties.particlesPerCell[0].z; 
	double e = 1.602176565e-19; // Elementary charge (Coulomb)  
	double K = 1e7 * e; // kinetic energy in Joule
 	double m = 1.672621777e-27; // Proton mass (kg) 
	double v_mod = universeProperties.speedOfLight / sqrt(1.0 + (m * universeProperties.speedOfLight * universeProperties.speedOfLight) / K);
	auto low = universeProperties.origin + 0.125 * universeSize;
	auto hig = low + 0.75 * universeSize;
	auto dist = distribution::uniform_pos_normal_speed<> ( 
			low, hig,
			Vector3<double> { 0, 0, 0 }, // mean value
			Vector3<double> { v_mod, v_mod, v_mod } // variance
	);
	auto universe = createUniverseFromDistribution(universeProperties, initProperties, numParticles, dist);

#ifdef ENABLE_DEBUG_OUTPUT
	// get the number of particles in all cells before the simulation begins for error checking
	//assert_decl(auto start_particles = countParticlesInDomain(universe));
#endif

	std::cout << "Running simulation..." << std::endl;

	// -- run the simulation --

	auto duration = simulateSteps(params.ncycles, universe);
	
	std::cout << "Simulation measurements: " << numParticles;
	std::cout << " initial particles, first step " << duration.firstStep << " seconds, " << (numParticles / duration.firstStep);
	std::cout << " pps, remaining steps " << duration.remainingSteps << " seconds, " << (numParticles*(params.ncycles-1))/duration.remainingSteps << " pps\n";

	// ----- finish ------

#ifdef ENABLE_DEBUG_OUTPUT
	// get the number of particles in all cells at the end of the simulation for error checking
	//assert_decl(auto end_particles = countParticlesInDomain(universe));
	//assert_eq(start_particles, end_particles) << "[Error]: Periodic boundary conditions on particles were not preserved!";
#endif

	std::cout << "Simulation finished successfully, producing output data..." << std::endl;

	std::string outputFilename = baseName + ".out";
	outputNumberOfParticlesPerCell(universe.cells, outputFilename);

	// be done
	return EXIT_SUCCESS;
}
