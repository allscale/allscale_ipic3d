#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
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

void traceParticle(Particle p, int T, const UniverseProperties& config, const InitProperties& init, int frame_interval, ParticleCount& res) {

	// extract some properties
	auto dt = config.dt;

	// get universe size
	auto universeSize = elementwiseProduct(config.cellWidth, config.size);
	auto universeLow  = config.origin;
	auto universeHigh = config.origin + universeSize;

	p.position += config.objectCenter;
	// exclude generated particles from the earth
	auto diff = p.position - config.objectCenter;
	double r2 = allscale::utils::sumOfSquares(diff);
	if (r2 < config.planetRadius * config.planetRadius) { 
		return;
	}

	for(std::size_t i=0; i<3; i++) {
		if (p.position[i] > universeHigh[i]) return;
		if (p.position[i] < universeLow[i])  return;
	}

	// simulate particle for the given number of time steps
	for(int t=0; t<T; t++) {

		// if the cell coordinate has changed
		auto pos = getCellCoordinates(config,p);

		// calculate 3 Cartesian components of the magnetic field
		double fac1 =  -init.externalMagneticField.z * pow(config.planetRadius, 3) / pow(allscale::utils::sumOfSquares(p.position), 2.5);
		Vector3<double> E, B;
		E = {0.0, 0.0, 0.0};
		B.x = 3.0 * p.position.x * p.position.z * fac1;
		B.y = 3.0 * p.position.y * p.position.z * fac1;
		B.z = (2.0 * pow(p.position.z, 2) - pow(p.position.x, 2) - pow(p.position.y, 2)) * fac1;
			
		// adaptive sub-cycling for computing velocity
		double B_mag = allscale::utils::sumOfSquares(B);
		double dt_sub = M_PI * config.speedOfLight / (4.0 * fabs(p.qom) * B_mag);
		int sub_cycles = int(dt / dt_sub) + 1;
		sub_cycles = std::min(sub_cycles, 100);
		dt_sub = dt / double(sub_cycles);

		for (int cyc_cnt = 0; cyc_cnt < sub_cycles; cyc_cnt++) {
			// update velocity
			p.updateVelocity(E,B,dt_sub);

			// update position
			p.updatePosition(dt_sub);
		}

		// support wrap-around -- periodic boundary conditions
		for(std::size_t i=0; i<3; i++) {
			if (p.position[i] > universeHigh[i]) p.position[i] -= universeSize[i];
			if (p.position[i] < universeLow[i])  p.position[i] += universeSize[i];
		}

		// remove particles from inside the sphere
		//   In fact, we just stop computing
		auto diff = p.position - config.objectCenter;
		double r2 = allscale::utils::sumOfSquares(diff);
		if (r2 <= config.planetRadius * config.planetRadius) {
			return;
		}

		// register particle in cell if necessary
		if (t % frame_interval == 0) {
			res.increment(t/frame_interval, pos);
		}
	}
}


int main(int argc, char** argv) {

	// ----- load and parse simulation parameters ------

	// parameters
	int N = 16*1000*1000;		// < number of particles
	int T = 150;		// < number of time steps
	int S = 10;		// < number of time steps between frames
	int R = 64;			// resolution of the result grid

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

	if (argc > 4) {
		R = atoi(argv[4]);
	}

	// some derived parameters
	int num_frames = T / S + 1;

	// print some introduction and summary information

	std::cout << "----- particle-in-cell tracer -----\n";
	std::cout << "Tracing " << N << " particles for " << T << " time steps in a " << R << "^3 grid recording a snapshot every " << S << " time steps ...\n";

	// set up relevant universe properties
	UniverseProperties config;
	config.dt = 0.01;
	config.speedOfLight = 299792458;
	config.size = { R, R, R };
	config.planetRadius = 6378137; // meter (Earth radius) 
	config.cellWidth = (20.0 / config.size.x) * config.planetRadius;
	config.FieldOutputCycle = 0;

	// these parameters are required for computations
	//config.objectCenter = {6.0 * config.planetRadius, 5.0 * config.planetRadius, 5.0 * config.planetRadius};
	config.objectCenter = {0.0, 0.0, 0.0};
	config.origin.x = config.objectCenter.x - config.size.x * config.cellWidth.x / 2.0; 
	config.origin.y = config.objectCenter.y - config.size.y * config.cellWidth.y / 2.0; 
	config.origin.z = config.objectCenter.z - config.size.z * config.cellWidth.z / 2.0; 

	InitProperties init;
	init.externalMagneticField = {0.0, 0.0, 3.07e-5};

	config.useCase = UseCase::Dipole;

	// run simulation
	auto start = std::chrono::high_resolution_clock::now();

	// the block size to be used
	int B = std::max(1000, N/1000);		// the blocking factor (for performance reasons)

	// values to initialize particles 
	double e = 1.602176565e-19; // Elementary charge (Coulomb)  
	double K = 1e7 * e; // kinetic energy in Joule
 	double m = 1.672621777e-27; // Proton mass (kg) 
	double v_mod = config.speedOfLight / sqrt(1.0 + (m * config.speedOfLight * config.speedOfLight) / K);

	// a map operator for a range of elements
	auto map = [=](int a, int b){
		ParticleCount res(num_frames,config.size);

		// create a random particle generator
		//double low0 = 1.25 * config.planetRadius, hig0 = 5.0 * config.planetRadius; 
		//auto R1 = Vector3<double>{low0, low0, low0};
		//auto R2 = Vector3<double>{hig0, hig0, hig0};
		auto low = config.origin + elementwiseProduct(config.cellWidth,config.size) / 8.0;
		auto hig = low + 0.75 * elementwiseProduct(config.cellWidth,config.size);
		distribution::uniform_pos_normal_speed<> next(
//				R1, R2, // within the universe
				low, hig,
				// speeds are normal distributed
				Vector3<double> {0.0, 0.0, 0.0},   // mean value
				Vector3<double> {v_mod, v_mod, v_mod}, // variance
				a*b
		);

		for(int i=a*B; i<b*B && i < N; i++) {
			// create a particle
			Particle p = next();
			p.q = e;
			p.qom = e / m;

			// trace its trajectory
			traceParticle(p,T+1,config,init,S,res);
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
	std::cout << "Throughput: " << ((T+1) * double(N)) / s << " particles/s \n";

	// save result
	{
		// use the current time to create a unique file name
		auto timeStamp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
		for(int t=0; t<num_frames; t++) {

			// generate file name
			char fileName[50];
			std::sprintf(fileName,"result_%ld.csv.%06d", timeStamp,t);

			// open file and dump results
			auto out = std::fstream(fileName, std::ios_base::out);
			out << "t,x,y,z,density\n";
			for(int x=0; x<config.size.x;x++) {
				for(int y=0; y<config.size.y;y++) {
					for(int z=0; z<config.size.z;z++) {
						double dx = x * config.cellWidth.x;
						double dy = y * config.cellWidth.y;
						double dz = z * config.cellWidth.z;
						out << t << "," << dx << "," << dy << "," << dz << "," << res.get(t,x,y,z) << "\n";
					}
				}
			}
		}
	}

	// be done
	return EXIT_SUCCESS;
}
