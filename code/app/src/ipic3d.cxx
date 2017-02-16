#include <cstdlib>
#include <iostream>

#include "ipic3d/app/cell.h"
#include "ipic3d/app/parameters.h"
#include "ipic3d/app/simulator.h"

#include "ipic3d/app/utils/points.h"

using namespace ipic3d;

Grid<Cell> initCells(const Parameters&);

Field initFields(const Parameters&);


int main(int argc, char** argv) {

	// ----- load and parse simulation parameters ------

	// check the passed arguments
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <config-file>\n";
		return EXIT_FAILURE;
	}

	// load input configuration
	std::cout << "Loading configuration file \"" << argv[1] << "\" ...\n";
	auto params = Parameters::read(argv[1]);


	params.dx = 10;
	params.dy = 10;
	params.dz = 10;

	// ----- initialize simulation environment ------

	// setup simulation
	std::cout << "Initializing simulation state ...\n";

	// initialize cells and contained particles
	auto cells = initCells(params);

	// initialize electric and magnetic force fields
	auto field = initFields(params);


	// ----- run the simulation ------

	// run simulation
	std::cout << "Running simulation ...\n";

	// extract size of grid
	auto size = cells.size();

	// get the time step
	double dt = params.dt;

	// print some infos for the user
	std::cout << "   Grid Size:       " << size[0] << " x " << size[1] << " x " << size[2] << "\n";
	std::cout << "   Time Step:       " << dt << "\n";
	std::cout << "   Number of steps: " << params.ncycles << "\n";

	// -- run the simulation --

	simulateSteps(params.ncycles, dt, cells, field);

	// ----- finish ------

	// be done
	std::cout << "Simulation finished successfully!\n";
	return EXIT_SUCCESS;
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
		cell.x = pos[0] * params.dx;
		cell.dx = params.dx;

		cell.y = pos[1] * params.dy;
		cell.dy = params.dy;

		cell.z = pos[2] * params.dz;
		cell.dz = params.dz;

		// -- add particles --

		// compute number of particles to be added
		int npcel = params.npcel[0];

		// add the requested number of parameters
		unsigned random_state = pos[0] * 10000 + pos[1] * 100 + pos[2];
		for (int i = 0; i < npcel; i++) {
			Particle p;

			p.x = cell.x + (rand_r(&random_state) / RAND_MAX) * cell.dx - cell.dx/2;
			p.y = cell.y + (rand_r(&random_state) / RAND_MAX) * cell.dy - cell.dy/2;
			p.z = cell.z + (rand_r(&random_state) / RAND_MAX) * cell.dz - cell.dz/2;

			// TODO: initialize the speed of particles
			p.dx = 0.0;
			p.dx = 0.0;
			p.dx = 0.0;

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

	// TODO: this is completely wrong => those values must not be shared!!
	pfor(zero, fieldSize,[&](const utils::Coordinate<3>& cur) {

		double x_displ, y_displ, z_displ, fac1;

		// init electrical field
		fields[cur].E.x = ebc[0];
		fields[cur].E.y = ebc[1];
		fields[cur].E.z = ebc[2];

		// init magnetic field
		fields[cur].B.x = params.B0x;
		fields[cur].B.y = params.B0y;
		fields[cur].B.z = params.B0z;

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
					fields[cur].Bext.x = 0.0;
					fields[cur].Bext.y = 0.0;
					fields[cur].Bext.z = 0.0;
				}

				break;
			}

			case UseCase::ParticleWave: {

				fields[cur].Bext.x = 0;
				fields[cur].Bext.y = 0;
				fields[cur].Bext.z = 0;

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

