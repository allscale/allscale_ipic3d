#include <cstdlib>
#include <iostream>

#include "allscale/api/user/data/vector.h"

#include "ipic3d/app/cell.h"
#include "ipic3d/app/parameters.h"
#include "ipic3d/app/simulator.h"
#include "ipic3d/app/universe.h"
#include "ipic3d/app/init_properties.h"

#include "ipic3d/app/utils/points.h"

using namespace ipic3d;

Grid<Cell> initCells(const InitProperties& initProperties, const UniverseProperties&);

Field initFields(const InitProperties& initProperties, const UniverseProperties&);

InitProperties initInitProperties(const Parameters& params);

Universe initUniverse(const Parameters&);

UniverseProperties initUniverseProperties(const Parameters&);

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
	auto universe = initUniverse(params);

	// ----- run the simulation ------

	std::cout << "Running simulation..." << std::endl;

	// -- run the simulation --

	simulateSteps(params.ncycles, universe);

	// ----- finish ------

	// be done
	std::cout << "Simulation finished successfully!" << std::endl;
	return EXIT_SUCCESS;
}

Universe initUniverse(const Parameters& params) {

	// initialize initial properties
	InitProperties initProperties = initInitProperties(params);

	std::cout << initProperties;

	// initialize universe properties
	UniverseProperties universeProperties = initUniverseProperties(params);

	std::cout << universeProperties;

	// create a universe with the given properties
	Universe universe(initCells(initProperties, universeProperties), initFields(initProperties, universeProperties), universeProperties);
	return universe;
}

InitProperties initInitProperties(const Parameters& params) {
	InitProperties properties;

	properties.numSteps = params.ncycles;

	for(int i = 0; i < params.ns; i++) {
		properties.driftVelocity.push_back({ params.u0[i], params.v0[i], params.w0[i] });
	}

	for(int i = 0; i < (params.ns+params.nstestpart); i++) {
		properties.particlesPerCell.push_back({ (unsigned)params.npcelx[i], (unsigned)params.npcely[i], (unsigned)params.npcelz[i] });
	}

	properties.magneticFieldAmplitude = {params.B0x, params.B0y, params.B0z};

	return properties;
}

UniverseProperties initUniverseProperties(const Parameters& params) {
	UniverseProperties properties;
	properties.cellWidth = { params.dx, params.dy, params.dz };
	properties.size = { params.nxc, params.nyc, params.nzc };
	properties.dt = params.dt;
	properties.useCase = params.useCase;
	return properties;
}

Grid<Cell> initCells(const InitProperties& initProperties, const UniverseProperties& properties) {

	const utils::Coordinate<3> zero = 0;							// a zero constant (coordinate [0,0,0])
	const utils::Coordinate<3> full = properties.size;				// a constant covering the full range


	// -- initialize the grid of cells --

	// the 3-D grid of cells
	Grid<Cell> cells(properties.size);								// the grid of cells containing the particles

	// -- initialize the state of each individual cell --

	// TODO: return this as a treeture
	allscale::api::user::pfor(zero, full, [&](const utils::Coordinate<3>& pos) {

		Cell& cell = cells[pos];

		// -- add particles --

		// compute number of particles to be added
		unsigned particlesPerCell = initProperties.particlesPerCell[0].x + initProperties.particlesPerCell[0].y + initProperties.particlesPerCell[0].z;

		// add the requested number of parameters
		unsigned random_state = pos[0] * 10000 + pos[1] * 100 + pos[2];
		for (unsigned i = 0; i < particlesPerCell; i++) {
			Particle p;

			Vector3<double> randVals = {(double)rand_r(&random_state) / RAND_MAX, (double)rand_r(&random_state) / RAND_MAX, (double)rand_r(&random_state) / RAND_MAX};
			// initialize particle position
			p.position = getCenterOfCell(pos, properties) + elementwiseProduct(randVals, properties.cellWidth) - properties.cellWidth/2;

			// TODO: initialize the speed of particles
			p.velocity = { 0, 0, 0 };

			// TODO: initialize charge
			p.q = 0.15;

			cell.particles.push_back(p);
		}

	});

	// return the initialized cells
	return std::move(cells);
}

Field initFields(const InitProperties& initProperties, const UniverseProperties& universeProperties) {

	using namespace allscale::api::user;

	// determine the field size
	utils::Size<3> zero = 0;
	utils::Size<3> fieldSize = universeProperties.size + coordinate_type(1);

	// the 3-D force fields
	Field fields(fieldSize);

	auto driftVel = initProperties.driftVelocity;
	assert_false(driftVel.empty()) << "Expected a drift velocity vector of at least length 1";
	auto ebc = crossProduct(driftVel[0], initProperties.magneticFieldAmplitude) * -1;

	pfor(zero, fieldSize,[&](const utils::Coordinate<3>& cur) {

		double fac1;

		// init electrical field
		fields[cur].E = ebc;

		// init magnetic field
		fields[cur].B = initProperties.magneticFieldAmplitude;

		// -- add earth model --

		switch(universeProperties.useCase) {

			case UseCase::Dipole: {

				// radius of the planet
				double a = universeProperties.objectRadius;

				auto objectCenter = universeProperties.objectCenter;
				auto location = getLocation(cur, universeProperties);

				auto diff = location - objectCenter;

				double r2 = allscale::api::user::data::sumOfSquares(diff);

				// Compute dipolar field B_ext

				if (r2 > a*a) {
					fac1 =  -universeProperties.magneticField.z * pow(a, 3) / pow(r2, 2.5);
					fields[cur].Bext.x = 3.0 * diff.x * diff.z * fac1;
					fields[cur].Bext.y = 3.0 * diff.y * diff.z * fac1;
					fields[cur].Bext.z = (2.0 * diff.z * diff.z - diff.x * diff.x - diff.y * diff.y) * fac1;
				} else { // no field inside the planet
					fields[cur].Bext = { 0.0, 0.0, 0.0 };
				}

				break;
			}

			case UseCase::ParticleWave: {

				fields[cur].Bext = { 0, 0, 0 };

				break;
			}

			default:
					assert_not_implemented()
						<< "The specified use case is not supported yet!";
		}

	});

	// return the produced field
	return std::move(fields);
}

