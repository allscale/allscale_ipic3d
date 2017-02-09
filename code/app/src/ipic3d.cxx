#include <cstdlib>
#include <iostream>

#include "allscale/api/user/data/grid.h"

#include "ipic3d/ipic3d.h"
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


SimulationState initializeSimulation(const Parameters&);

void runSimulation(SimulationState&);


int main(int argc, char** argv) {

	// check the passed arguments
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <config-file>\n";
		return EXIT_FAILURE;
	}

	// load input configuration
	auto params = Parameters::read(argv[1]);


	// setup simulation
	auto state = initializeSimulation(params);

	// run simulation
	runSimulation(state);

	// be done
	return EXIT_SUCCESS;
}


// initialize grid of cells
// TODO: is this the right place for this routine?
void initGrid(const utils::Coordinate<3>& zero, const utils::Coordinate<3>& full, allscale::api::user::data::Grid<Cell,3> &cells, const Parameters& params){
	// TODO: revise
	allscale::api::user::pfor(zero, full, [&](const utils::Coordinate<3>& pos) {

		Cell& cell = cells[pos];

		// physical properties
		cell.x = pos[0] * params.dx;
		cell.dx = params.dx;

		cell.y = pos[1] * params.dy;
		cell.dy = params.dy;

		cell.z = pos[2] * params.dz;
		cell.dz = params.dz;

	});
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



SimulationState initializeSimulation(const Parameters& params) {

	// ----- setup -----
	utils::Size<3> size = {params.nxc, params.nyc, params.nzc};		// the size of the grid

	const utils::Coordinate<3> zero = 0;							// a zero constant (coordinate [0,0,0])
	const utils::Coordinate<3> full = size;							// a constant covering the full range


	// the 3-D grid of cells
	Grid<Cell,3> cells(size);									// the grid of cells containing the particles

	// initialize all cells
	initGrid(zero, full, cells, params);

	// the field size
	// 		fields are defined on nodes
	utils::Size<3> fieldSize = {params.nxc + 1, params.nyc + 1, params.nzc + 1};

	// the 3-D force field
	Field field(fieldSize);

	// TODO: initFields();
	initDipole(field, params);

	// the 3-D density field
	// rho is charge density
	// J is current density
	// TODO:		we may need to compute them twice: on centers of cells and on nodes
	Density density(fieldSize);
	// initDensity();

	// a 3-D structure collecting contributions of cells to the density
	// TODO: revise, e.g. why do we need fieldSize * 2 DensityCells
	Grid<DensityCell,3> densityContributions(fieldSize*2);

	// initialize contributions to 0
	allscale::api::user::pfor(fieldSize*2,[&](const utils::Coordinate<3>& pos){
		auto& entry = densityContributions[pos];
		entry.rho = 0.0;
		entry.J = 0.0;
	});

	// insert particles in all cells
	int npcel = params.npcel[0];
	allscale::api::user::pfor(zero, full, [&](const utils::Coordinate<3>& pos) {

		Cell& cell = cells[pos];

		// physical properties
		for (int i = 0; i < npcel; i++) {
			double x, y, z;
			x = cell.x + (rand() / RAND_MAX) * cell.dx;
			y = cell.y + (rand() / RAND_MAX) * cell.dy;
			z = cell.z + (rand() / RAND_MAX) * cell.dz;
			Particle p{x, y, z, 0.0, 0.0, 0.0, 0.0};
			cell.particles.push_back(p);
		}
	});

	// initialize all particles
	//double dt = 0.01, tcur = 0.0, tend = 230.0;
	double dt = 0.01, tcur = 0.0, tend = 6.0;
	allscale::api::user::pfor(zero, size, [&](const utils::Coordinate<3>& pos) {
		cells[pos].initParticles(pos, dt, true);
	});

//	/*// create the output file
//	auto& manager = utils::io::FileIOManager::getInstance();
//	// create the result file
//	// text file
//	auto logFile = manager.createEntry("ipic3d_particles_dipole_log.txt", utils::io::Mode::Text);
//	auto& outtxt = manager.openOutputStream(logFile);
//	outtxt << "x,y,z,vx,vy,vz\n";*/

	// return the created simulation state
	return { std::move(cells), std::move(field) };
}

void runSimulation(SimulationState& state) {

	auto& cells = state.cells;
	auto& field = state.field;

	auto size = cells.size();
	const utils::Coordinate<3> zero = 0;							// a zero constant (coordinate [0,0,0])

	Grid<std::vector<Particle>,3> particleTransfers(size * 3);	// a grid of buffers for transferring particles between cells

	// initialize all particles
	double dt = 0.00625, tcur = 0.0, tend = 4.0;
	allscale::api::user::pfor(zero, size, [&](const utils::Coordinate<3>& pos) {
		cells[pos].initParticles(pos, dt, false);
	});


	// propagate particles
	while( tcur <= tend ) {

		// TODO: interpolation of particles to grid: project particles to the density field
		/*allscale::api::user::pfor(zero,full,[&](const utils::Coordinate<3>& pos) {
			cells[pos].interP2G(pos, densityContributions);
		});

		// TODO: update density field
		allscale::api::user::pfor(zero,fieldSize,[&](const utils::Coordinate<3>& pos) {
			DensityCell& entry = density[pos];

			// reset density field
			entry.rho = 0.0;
			entry.J = 0.0;

			// aggregate contributions of adjacent cell
			for(int i=0; i<2; i++) {
				for(int j=0; j<2; j++) {
					for(int k=0; k<2; k++) {
						auto& cur = densityContributions[(pos * 2) + utils::Coordinate<3>{i,j,k}];
						entry.rho += cur.rho;
						entry.J += cur.J;
					}
				}
			}
		});


		// TODO: propagate field forces according to Maxwell's equations (GMRES)
		FieldSolver();*/

		// move particles
		allscale::api::user::pfor(zero,size,[&](const utils::Coordinate<3>& pos){
			cells[pos].BorisMover(pos, field, particleTransfers, dt, true);
		});
		// -- implicit global sync - TODO: can this be eliminated? --

		// increment time step
		tcur += dt;
	}

}
