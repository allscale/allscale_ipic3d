#pragma once

#include <chrono>

#include "ipic3d/app/vector.h"
#include "ipic3d/app/universe.h"
#include "ipic3d/app/simulator.h"

#include "ipic3d/app/mpi_context.h"

namespace ipic3d {

	using namespace allscale::utils;

	const uint64_t NUM_TIME_STEPS = 100;
	const double DELTA_T = 0.15;

	const Vector3<int64_t> GRID_SIZE{ 64, 64, 64 };

	const Vector3<double> CELL_WIDTH{ 10, 10, 10 };

	const Vector3<double> UNIVERSE_SIZE{
		GRID_SIZE.x * CELL_WIDTH.x,
		GRID_SIZE.y * CELL_WIDTH.y,
		GRID_SIZE.z * CELL_WIDTH.z
	};

	int processUniverse(Universe& universe, std::uint64_t numParticles, std::uint64_t numTimeSteps) {
		const bool dumpParticles = std::getenv("DUMP_PARTICLE_POSITION");

		// --- sample output initial state --
		if(dumpParticles) outputParticlePositions(universe.cells, std::string("t_begin.txt"));


		// ----- run the simulation ------
		if (MPI_Context::isMaster()) std::cout << "Running simulation..." << std::endl;

		auto start = std::chrono::high_resolution_clock::now();
		simulateSteps(numTimeSteps, universe);
		auto end = std::chrono::high_resolution_clock::now();

		if (MPI_Context::isMaster()) std::cout << "Simulation Finished" << std::endl;

		double s = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
		if (MPI_Context::isMaster()) std::cout << "Simulation took " << s << "s\n";
		if (MPI_Context::isMaster()) std::cout << "Throughput: " << (numTimeSteps * numParticles) / s << " particles/s \n";

		if(dumpParticles) outputParticlePositions(universe.cells, std::string("t_end.txt"));

		return EXIT_SUCCESS;
	}

	template<typename Distribution>
	int processDistribution(const Distribution& dist, std::uint64_t numParticles, std::uint64_t numTimeSteps) {

		// set up universe
		UniverseProperties universeProperties;
		universeProperties.cellWidth = CELL_WIDTH;
		universeProperties.size = GRID_SIZE;
		universeProperties.dt = DELTA_T;
		universeProperties.FieldOutputCycle = 0;

		// set up init properties
		InitProperties initProperties;
		initProperties.driftVelocity.push_back(0);

		// initialize universe
		if (MPI_Context::isMaster()) std::cout << "Creating Particles ..." << std::endl;
		auto universe = createUniverseFromDistribution(universeProperties, initProperties, numParticles, dist);

		// -- run the simulation --

		return processUniverse(universe, numParticles, numTimeSteps);

	}

	int processUniform(std::uint64_t numParticles, std::uint64_t numTimeSteps) {
		return processDistribution(
			ipic3d::distribution::uniform<>(
				Vector3<double> { 0, 0, 0 },
				UNIVERSE_SIZE,
				Vector3<double> { -0.2, -0.2, -0.2},
				Vector3<double> { +0.2, +0.2, +0.2}
		),
			numParticles,
			numTimeSteps
			);
	}

	int processCluster(std::uint64_t numParticles, std::uint64_t numTimeSteps) {
		return processDistribution(
			ipic3d::distribution::normal<>(
				UNIVERSE_SIZE / 2,
				UNIVERSE_SIZE / 5,
				Vector3<double> { -0.2, -0.2, -0.2},
				Vector3<double> { +0.2, +0.2, +0.2}
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
				UNIVERSE_SIZE / 2,
				UNIVERSE_SIZE.x / 10
			),
			numParticles,
			numTimeSteps
		);
	}

	int processBeam(std::uint64_t numParticles, std::uint64_t numTimeSteps) {
		return processDistribution(
			ipic3d::distribution::normal<>(
				UNIVERSE_SIZE / 100,
				UNIVERSE_SIZE / 500,
				UNIVERSE_SIZE / numTimeSteps * 0.95 / DELTA_T,
				UNIVERSE_SIZE / numTimeSteps * 0.95 / DELTA_T
				),
			numParticles,
			numTimeSteps
		);
	}

	int processBenchmark(const std::string& config) {

		if (MPI_Context::isMaster()) std::cout << "Processing benchmark " << config << "\n";

		// utility for help messages
		auto showHelp = [&]() {
			if (MPI_Context::isMaster()) std::cout << "Benchmark designation: :X:N\n"
				<< "      X ... benchmark type:\n"
				<< "           U ... uniform\n"
				<< "           C ... cluster\n"
				<< "           E ... explosion\n"
				<< "           B ... beam\n"
				<< "      N ... total number of particles\n";
			return EXIT_FAILURE;
		};

		// check the format
		if(config.length() <= 3 || config[0] != ':' || config[2] != ':') {
			return showHelp();
		}

		// all benchmarks comprise 100 time steps
		std::uint64_t numTimeSteps = NUM_TIME_STEPS;

		// parse number of particles
		std::uint64_t numParticles = std::atoll(config.substr(3, config.length()).c_str());

		// dispatch based on use case
		switch(config[1]) {
			case 'U': return processUniform(numParticles, numTimeSteps);
			case 'C': return processCluster(numParticles, numTimeSteps);
			case 'E': return processExplosion(numParticles, numTimeSteps);
			case 'B': return processBeam(numParticles, numTimeSteps);
			default:  return showHelp();
		}
	}

} // end namespace ipic3d
