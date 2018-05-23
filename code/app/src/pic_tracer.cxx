#include <chrono>
#include <cstdlib>
#include <iostream>


#include "allscale/api/user/algorithm/pfor.h"

#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/particle.h"
#include "ipic3d/app/universe_properties.h"

using namespace allscale::utils;
using namespace allscale::api::user::algorithm;
using namespace ipic3d;


void traceParticle(Particle p, int T, const UniverseProperties& config, const Field& field) {

	// extract some properties
	auto dt = config.dt;
	auto cellWidth = config.cellWidth;

	double vol = config.cellWidth.x * config.cellWidth.y * config.cellWidth.z;

	// get universe size
	auto universeSize = elementwiseProduct(config.cellWidth, config.size);
	auto universeLow  = config.origin;
	auto universeHigh = config.origin + universeSize;

	// get the currently active cell position
	utils::Coordinate<3> pos = { -1, -1, -1 };
	auto cellOrigin = getOriginOfCell(pos, config);

	// initialize extracted field parameters
	Vector3<double> Es[2][2][2];
	Vector3<double> Bs[2][2][2];

	// simulate particle for the given number of time steps
	for(int t=0; t<T; t++) {

		// if the cell coordinate has changed
		auto curPos = getCellCoordinates(config,p);
		if (pos != curPos) {

			// update position and cell origin
			pos = curPos;
			cellOrigin = getOriginOfCell(pos,config);

			// update Es and Bs based on current cell position
			for(int i=0; i<2; i++) {
				for(int j=0; j<2; j++) {
					for(int k=0; k<2; k++) {
						utils::Coordinate<3> cur({pos[0]+i+1,pos[1]+j+1,pos[2]+k+1});
						Es[i][j][k] = field[cur].E;
						Bs[i][j][k] = field[cur].B;
					}
				}
			}
		}

		// get the fractional distance of the particle from the cell origin
		const auto relPos = allscale::utils::elementwiseDivision((p.position - cellOrigin), cellWidth);

		// interpolate forces
		auto E = trilinearInterpolationF2P(Es, relPos, vol);
		auto B = trilinearInterpolationF2P(Bs, relPos, vol);

		// update velocity
		p.updateVelocity(E,B,dt);

		// update position
		p.updatePosition(dt);

		// support wrap-around
		for(std::size_t i=0; i<3; i++) {
			if (p.position[i] > universeHigh[i]) p.position[i] -= universeSize[i];
			if (p.position[i] < universeLow[i])  p.position[i] += universeSize[i];
		}
	}
}


int main(int argc, char** argv) {

	// ----- load and parse simulation parameters ------

	// parameters
	int N = 1000 * 1000;		// < number of particles
	int T = 1000;		// < number of time steps

	// take command line parameters
	if (argc > 1) {
		N = atoi(argv[1]);
	}

	if (argc > 2) {
		T = atoi(argv[2]);
	}

	// print some introduction and summary information

	std::cout << "----- particle-in-cell tracer -----\n";
	std::cout << "Tracing " << N << " particles for " << T << " time steps ...\n";

	// set up relevant universe properties
	UniverseProperties config;
	config.dt = 0.01;
	config.cellWidth = 10;
	config.size = { 64, 64, 64 };
	config.FieldOutputCycle = 0;

	config.useCase = UseCase::Dipole;

	// create a field
	InitProperties initProps;
	initProps.driftVelocity.push_back(0);
	Field field = initFields(initProps, config);

	// run simulation
	auto start = std::chrono::high_resolution_clock::now();
	pfor(0,N,[&,T](int){

		// create a particle
		Particle p;		// TODO: define initial position and velocity
		p.position = {  1,  1,  1 };
		p.velocity = { 10, 10, 10 };

		// trace its trajectory
		traceParticle(p,T,config,field);
	});
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Simulation Finished" << std::endl;

	// done
	double s = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / (1e9);
	std::cout << "Simulation took " << s << "s\n";
	std::cout << "Throughput: " << (T * N) / s << " particles/s \n";

	// be done
	return EXIT_SUCCESS;
}
