#include <cstdlib>
#include <iostream>

#include "ipic3d/app/cell.h"
#include "ipic3d/app/parameters.h"
#include "ipic3d/app/simulator.h"
#include "ipic3d/app/universe.h"

#include "ipic3d/app/utils/points.h"

using namespace ipic3d;

Grid<Cell> initCells(const Parameters&);

Field initFields(const Parameters&);

Universe initUniverse(const Parameters&);

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

	// run simulation
	std::cout << "Running simulation ..." << std::endl;

	// extract size of grid
	auto size = universe.size();

	// get the time step
	double dt = params.dt;

	// print some infos for the user
	std::cout << "   Grid Size:       " << size[0] << " x " << size[1] << " x " << size[2] << std::endl;
	std::cout << "   Time Step:       " << dt << std::endl;
	std::cout << "   Number of steps: " << params.ncycles << std::endl;

	// -- run the simulation --

	simulateSteps(params.useCase, params.ncycles, dt, universe);

	// ----- finish ------

	// be done
	std::cout << "Simulation finished successfully!" << std::endl;
	return EXIT_SUCCESS;
}

Universe initUniverse(const Parameters& params) {
	return Universe(initCells(params), initFields(params));
}

Grid<Cell> initCells(const Parameters& params) {

	// ----- setup -----
	utils::Size<3> size = {params.nxc, params.nyc, params.nzc};		// the size of the grid

	const utils::Coordinate<3> zero = 0;							// a zero constant (coordinate [0,0,0])
	const utils::Coordinate<3> full = size;							// a constant covering the full range


	// -- initialize the grid of cells --

	// the 3-D grid of cells
	Grid<Cell> cells(size);									// the grid of cells containing the particles

	// -- initialize the state of each individual cell --

	// TODO: return this as a treeture
	allscale::api::user::pfor(zero, full, [&](const utils::Coordinate<3>& pos) {

		Cell& cell = cells[pos];

		// physical properties
		cell.spacing = {params.dx, params.dy, params.dz};
		Vector3<double> tempPos = { (double)pos[0], (double)pos[1], (double)pos[2] };
		cell.center = elementwiseProduct(tempPos, cell.spacing);

		// -- add particles --

		// compute number of particles to be added
		int npcel = params.npcel[0];

		// add the requested number of parameters
		unsigned random_state = pos[0] * 10000 + pos[1] * 100 + pos[2];
		for (int i = 0; i < npcel; i++) {
			Particle p;

			Vector3<double> randVals = {(double)rand_r(&random_state) / RAND_MAX, (double)rand_r(&random_state) / RAND_MAX, (double)rand_r(&random_state) / RAND_MAX};
			p.position = cell.center + elementwiseProduct(randVals, cell.spacing) - cell.spacing / 2.0;

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

Field initFields(const Parameters& params) {

	using namespace allscale::api::user;

	// determine the field size
	utils::Size<3> zero = 0;
	utils::Size<3> fieldSize = {params.nxc + 1, params.nyc + 1, params.nzc + 1};

	// the 3-D force fields
	Field fields(fieldSize);

	// TODO: use one kind of vector everywhere
	data::Vector<double,3> a{params.u0[0], params.v0[0], params.w0[0]};
	data::Vector<double,3> b{params.B0x, params.B0y, params.B0z};
	auto ebc = crossProduct(a,b) * -1;

	pfor(zero, fieldSize,[&](const utils::Coordinate<3>& cur) {

		double x_displ, y_displ, z_displ, fac1;

		// init electrical field
		fields[cur].E = { ebc[0], ebc[1], ebc[2] };

		// init magnetic field
		fields[cur].B = { params.B0x, params.B0y, params.B0z };

		// -- add earth model --

		switch(params.useCase) {

			case UseCase::Dipole: {

				// radius of the planet
				double a = params.L_square;

				double xc = params.x_center;
				double yc = params.y_center;
				double zc = params.z_center;

				double x = cur[0];
				double y = cur[1];
				double z = cur[2];

				double r2 = ((x-xc)*(x-xc)) + ((y-yc)*(y-yc)) + ((z-zc)*(z-zc));

				// Compute dipolar field B_ext

				if (r2 > a*a) {
					x_displ = x - xc;
					y_displ = y - yc;
					z_displ = z - zc;
					fac1 =  -params.B1z * a * a * a / pow(r2, 2.5);
					fields[cur].Bext.x = 3.0 * x_displ * z_displ * fac1;
					fields[cur].Bext.y = 3.0 * y_displ * z_displ * fac1;
					fields[cur].Bext.z = (2.0 * z_displ * z_displ - x_displ * x_displ - y_displ * y_displ) * fac1;
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

