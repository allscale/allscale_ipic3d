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

///**
// * This function updates the position of all particles within a cell for a single
// * time step, considering the given field as a driving force.
// *
// * @param universeProperties the properties of this universe
// * @param cell the cell whose particles are moved
// * @param pos the coordinates of this cell in the grid
// * @param field the most recently computed state of the surrounding force fields
// */
//void moveParticles(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field) {
//
//	assert_true(pos.dominatedBy(universeProperties.size)) << "Position " << pos << " is outside universe of size " << universeProperties.size;
//
//	// quick-check
//	if (cell.particles.empty()) return;
//
//	// -- move the particles in space --
//
//	// extract forces
//	// TODO: move this to some C++ structure
//	Vector3<double> Es[2][2][2];
//	Vector3<double> Bs[2][2][2];
//	for(int i=0; i<2; i++) {
//		for(int j=0; j<2; j++) {
//			for(int k=0; k<2; k++) {
//				utils::Coordinate<3> cur({pos[0]+i+1,pos[1]+j+1,pos[2]+k+1});
//				Es[i][j][k] = field[cur].E;
//				Bs[i][j][k] = field[cur].B;
//			}
//		}
//	}
//
//	const auto cellOrigin = getOriginOfCell(pos, universeProperties);
//
//	double vol = universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;
//
//	// update particles
////		allscale::api::user::algorithm::pfor(cell.particles, [&](Particle& p){
//	for(std::size_t i=0; i<cell.particles.size(); ++i) {
//		Particle& p = cell.particles[i];
//		// Docu: https://www.particleincell.com/2011/vxb-rotation/
//		// Code: https://www.particleincell.com/wp-content/uploads/2011/07/ParticleIntegrator.java
//
//		// get the fractional distance of the particle from the cell origin
//		const auto relPos = allscale::utils::elementwiseDivision((p.position - cellOrigin), (universeProperties.cellWidth));
//
//		// interpolate
//		auto E = trilinearInterpolationF2P(Es, relPos, vol);
//		auto B = trilinearInterpolationF2P(Bs, relPos, vol);
//
//		// update velocity
//		p.updateVelocity(E, B, universeProperties.dt);
//
//		// update position
//		p.updatePosition(universeProperties.dt);
//	}
////		});
//
//}



void traceParticle(int T, const UniverseProperties& config, const Field& field) {

	// extract some properties
	auto dt = config.dt;
	auto cellWidth = config.cellWidth;

	// create a particle
	Particle p;		// TODO: define initial position and velocity

	double vol = config.cellWidth.x * config.cellWidth.y * config.cellWidth.z;

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

	// create a field
	InitProperties initProps;
	initProps.driftVelocity.push_back(0);
	Field field = initFields(initProps, config);

	// run simulation
	auto start = std::chrono::high_resolution_clock::now();
	pfor(0,N,[&,T](int){
		// simulate one particle
		traceParticle(T,config,field);
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
