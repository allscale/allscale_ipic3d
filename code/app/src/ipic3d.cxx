#include <cstdlib>
#include <iostream>

#include "allscale/api/user/data/grid.h"

#include "ipic3d/cell.h"
#include "ipic3d/parameters.h"
#include "ipic3d/utils/points.h"

using namespace ipic3d;
using namespace allscale::api::user::data;



struct SimulationState {

	// the current state of the cells
	Grid<Cell,3> cells;

	// the current state of the field
	Field field;

};


Grid<Cell,3> initCells(const Parameters&);

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



	// ----- initialize simulation environment ------

	// setup simulation
	std::cout << "Initializing simulation state ...\n";
	auto cells = initCells(params);
	auto field = initFields(params);



	// ----- run the simulation ------

	// run simulation
	std::cout << "Running simulation ...\n";

	// TODO: those values for sure need to be extracted from the parameters
	double dt = params.dt, tcur = 0.0, tend = 4.0;
	std::cout << "   Start Time: " << tcur << "\n";
	std::cout << "   End Time:   " << tend << "\n";
	std::cout << "   Time Step:  " << dt << "\n";

	// extract size of grid
	const utils::Coordinate<3> zero = 0;						// a zero constant (coordinate [0,0,0])
	auto size = cells.size();

	Grid<std::vector<Particle>,3> particleTransfers(size * 3);	// a grid of buffers for transferring particles between cells

	// initialize all particles
	allscale::api::user::pfor(zero, size, [&](const utils::Coordinate<3>& pos) {
		cells[pos].initParticles(pos, dt, false);
	});

	// run time loop for the simulation
	while( tcur <= tend ) {

		// move particles
		allscale::api::user::pfor(zero,size,[&](const utils::Coordinate<3>& pos){
			cells[pos].BorisMover(pos, field, particleTransfers, dt, true);
		});
		// -- implicit global sync - TODO: can this be eliminated? --

		// increment time step
		tcur += dt;
	}


	// ----- finish ------

	// be done
	std::cout << "Simulation finished successfully!\n";
	return EXIT_SUCCESS;
}


// initialize fields for the 3D DIPOLE simulation
// TODO: is this the right place for this routine?
void initDipole(Field &fields, const Parameters& params){

	allscale::api::user::data::Vector<double,3> a{params.u0[0], params.v0[0], params.w0[0]};
	allscale::api::user::data::Vector<double,3> b{params.B0x, params.B0y, params.B0z};
	auto ebc = crossProduct(a,b) * -1;

	double x_displ, y_displ, z_displ, fac1;
	allscale::api::user::pfor(fields.size(),[&](const utils::Coordinate<3>& cur) {
		fields[cur].E.x = ebc[0];
		fields[cur].E.y = ebc[1];
		fields[cur].E.z = ebc[2];

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
		fields[cur].B.x = params.B0x;// + Bx_ext[i][j][k]
		fields[cur].B.y = params.B0y;// + By_ext[i][j][k]
		fields[cur].B.z = params.B0z;// + Bz_ext[i][j][k]
	});

	// TODO: we need to compute Bc on centers of each node
	// 		that means we need to aggregate all the eight values from the nodes of each cell

	// TODO: communicateCenterBC_P

	// TODO: compute rho also with interpN2C

}


Grid<Cell,3> initCells(const Parameters& params) {

	// ----- setup -----
	utils::Size<3> size = {params.nxc, params.nyc, params.nzc};		// the size of the grid

	const utils::Coordinate<3> zero = 0;							// a zero constant (coordinate [0,0,0])
	const utils::Coordinate<3> full = size;							// a constant covering the full range


	// -- initialize the grid of cells --

	// the 3-D grid of cells
	Grid<Cell,3> cells(size);									// the grid of cells containing the particles

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

		// add particles
		int npcel = params.npcel[0];
		for (int i = 0; i < npcel; i++) {
			double x, y, z;
			x = cell.x + (rand() / RAND_MAX) * cell.dx;
			y = cell.y + (rand() / RAND_MAX) * cell.dy;
			z = cell.z + (rand() / RAND_MAX) * cell.dz;
			Particle p{x, y, z, 0.0, 0.0, 0.0, 0.0};
			cell.particles.push_back(p);
		}

		// initialize all the particles
		cell.initParticles(pos, params.dt, true);
	});

	// return the initialized cells
	return std::move(cells);
}

Field initFields(const Parameters& params) {

	using namespace allscale::api::user;

	// determine the field size
	utils::Size<3> fieldSize = {params.nxc + 1, params.nyc + 1, params.nzc + 1};

	// the 3-D force fields
	Field fields(fieldSize);


	data::Vector<double,3> a{params.u0[0], params.v0[0], params.w0[0]};
	data::Vector<double,3> b{params.B0x, params.B0y, params.B0z};
	auto ebc = crossProduct(a,b) * -1;

	double x_displ, y_displ, z_displ, fac1;
	pfor(fieldSize,[&](const utils::Coordinate<3>& cur) {
		fields[cur].E.x = ebc[0];
		fields[cur].E.y = ebc[1];
		fields[cur].E.z = ebc[2];

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
		fields[cur].B.x = params.B0x;// + Bx_ext[i][j][k]
		fields[cur].B.y = params.B0y;// + By_ext[i][j][k]
		fields[cur].B.z = params.B0z;// + Bz_ext[i][j][k]
	});

	// TODO: we need to compute Bc on centers of each node
	// 		that means we need to aggregate all the eight values from the nodes of each cell

	// TODO: communicateCenterBC_P

	// TODO: compute rho also with interpN2C


	// TODO: initFields();

	// return the produced field
	return std::move(fields);
}
