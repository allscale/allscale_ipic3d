#include <chrono>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/user/algorithm/preduce.h"

#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/particle.h"
#include "ipic3d/app/universe_properties.h"

using namespace allscale::utils;
using namespace allscale::api::user::algorithm;
using namespace ipic3d;

// --- density tracing ---

/**
 * The class aggregating the resulting particle density for one time step.
 */
class ParticleCount {

	int num_frames;
	utils::Coordinate<3> size;
	std::vector<std::uint32_t> counts;

public:

	ParticleCount() {}

	ParticleCount(const ParticleCount&) = delete;
	ParticleCount(ParticleCount&&) = default;

	ParticleCount(int num_frames, const utils::Coordinate<3>& size)
		: size(size), counts(num_frames*size.x*size.y*size.z,0) {}

	ParticleCount& operator=(const ParticleCount&) = delete;
	ParticleCount& operator=(ParticleCount&&) = default;

	std::uint32_t get(int t, int x, int y, int z) const {
		return counts[flatten(t,x,y,z)];
	}

	void increment(int t, int x, int y, int z, int num = 1) {
		counts[flatten(t,x,y,z)] += num;
	}

	void increment(int t, const utils::Coordinate<3>& pos, int num = 1) {
		increment(t,pos.x,pos.y,pos.z,num);
	}

	ParticleCount& operator+=(const ParticleCount& other) {
		assert_eq(size,other.size);
		for(std::size_t i=0; i<counts.size(); i++) {
			counts[i] += other.counts[i];
		}
		return *this;
	}

private:

	std::uint32_t flatten(int t, int x, int y, int z) const {
		return ((t*size.x + x)*size.y+y)*size.z+z;
	}

};



// --- tracing particles ---


void traceParticle(Particle p, int T, const UniverseProperties& config, const InitProperties& fieldProperties, int frame_interval, ParticleCount& res) {

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
						Vector3<double> loc = universeLow + elementwiseProduct(config.cellWidth, cur);
						auto fieldNode = getDipoleFieldAt(loc, fieldProperties, config);
						Es[i][j][k] = fieldNode.E;
						Bs[i][j][k] = fieldNode.B;
					}
				}
			}
		}

		// register particle in cell if necessary
		if (t % frame_interval == 0) {
			res.increment(t/frame_interval, pos);
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
	int T = 1000;				// < number of time steps
	int S = 1000;				// < number of time steps between frames

	// take command line parameters
	if (argc > 1) {
		N = atoi(argv[1]);
	}

	if (argc > 2) {
		T = atoi(argv[2]);
	}

	if (argc > 3) {
		S = atoi(argv[3]);
	}

	// some derived parameters
	int num_frames = T / S;

	// print some introduction and summary information

	std::cout << "----- particle-in-cell tracer -----\n";
	std::cout << "Tracing " << N << " particles for " << T << " time steps recording a snapshot every " << S << " time steps ...\n";

	// set up relevant universe properties
	UniverseProperties config;
	config.dt = 0.01;
	config.cellWidth = 10;
//	config.size = { 64, 64, 64 };
	config.size = { 16, 16, 16 };
	config.FieldOutputCycle = 0;

	config.useCase = UseCase::Dipole;
	config.magneticField = { 0, 0, 1 };

	// create a field
	InitProperties initProps;
	initProps.driftVelocity.push_back(0);

	// run simulation
	auto start = std::chrono::high_resolution_clock::now();

	// Old pfor-based version in case AllScale compiler will be unable to compile reduction
//		pfor(0,N,[&,T,config,initProps](int){
//
//			// create a particle
//			Particle p;		// TODO: define initial position and velocity
//			p.position = {  1,  1,  1 };
//			p.velocity = { 10, 10, 10 };
//
//			// trace its trajectory
//			traceParticle(p,T,config,initProps);
//		});


	// the block size to be used
	int B = 1000;		// the blocking factor ..

	// a map operator for a range of elements
	auto map = [=](int a, int b){
		ParticleCount res(num_frames,config.size);

		// create a random particle generator
		auto low = config.origin;
		auto hig = low + elementwiseProduct(config.cellWidth,config.size);
		distribution::uniform<> next(
				low,hig, // within the universe
				// speeds are constant
				Vector3<double> { -0.2, -0.2, -0.2},
				Vector3<double> { +0.2, +0.2, +0.2},
				a*b
		);

		for(int i=a*B; i<b*B && i < N; i++) {
			// create a particle
			Particle p = next();

			// trace its trajectory
			traceParticle(p,T,config,initProps,S,res);
		}
		return res;
	};

	// the reduction operation
	auto reduce = [](ParticleCount&& a, const ParticleCount& b) {
		a += b;
		return std::move(a);
	};

	auto res = preduce(0,N/B+1,map,reduce).get();

	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Simulation Finished" << std::endl;

	// print performance summary
	double s = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / (1e9);
	std::cout << "Simulation took " << s << "s\n";
	std::cout << "Throughput: " << (T * N) / s << " particles/s \n";

	// save result
	for(int t=0; t<num_frames; t++) {
		std::cout << "t=" << t << "\n";
		for(int x=0; x<config.size.x;x++) {
			for(int y=0; y<config.size.y;y++) {
				for(int z=0; z<config.size.z;z++) {
					std::cout << res.get(0,x,y,z) << " ";
				}
				std::cout << "\n";
			}
			std::cout << "\n";
		}
	}

	// be done
	return EXIT_SUCCESS;
}
