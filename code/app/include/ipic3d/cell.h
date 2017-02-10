#pragma once

#include <vector>

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/operator/pfor.h"

#include "ipic3d/particle.h"
#include "ipic3d/field_node.h"
#include "ipic3d/utils/points.h"
#include "ipic3d/utils/static_grid.h"

namespace ipic3d {

	// Earth parameters
	static const double Re = 6378137.0; 		// meter (Earth radius)
	static const double B0 = 3.07e-5; 			// Tesla
	// Other parameters
	static const double e = 1.602176565e-19; 	// Elementary charge (Coulomb)
	static const double m = 1.672621777e-27; 	// Proton mass (kg)
	static const double c = 299792458.0; 		// speed of light (m/s)
	static const double K = 1e7 * e;    		// kinetic energy in eV converted to Joules


	struct DensityCell {
		double rho;			// charge density
		double J;			// current density
	};

	using Density = allscale::api::user::data::Grid<DensityCell,3>;	// a 3D grid of density cells

	/**
	 * The structure of a single cell, forming a container for a set of particles
	 * located within a confined area of space.
	 */
	struct Cell {

		using Coord = utils::Coordinate<3>;

		// the cell position
		double x, y, z;

		// the cell grid spacing
		double dx, dy, dz;

		// the list of local particles
		std::vector<Particle> particles;

		/**
		 * Requires this cell, being located at position `pos`, to project the effect
		 * of its contained particles to the density grid. Contributions are stored
		 * within the given contributions grid
		 */
		void projectToDensityField(const Coord& pos, allscale::api::user::data::Grid<DensityCell,3>& contributions) const {

			// quick-check
			if (particles.empty()) return;		// nothing to contribute

			// init aggregated densities of neighboring cells
			DensityCell res[2][2][2];
			for(int i=0; i<2; i++) {
				for(int j=0; j<2; j++) {
					for(int k=0; k<2; k++) {
						res[i][j][k].rho = 0.0;
						res[i][j][k].J = 0.0;
					}
				}
			}

			// aggregate particles
			// TODO data race on res, should be avoided
			for(const auto& p : particles) {
				for(int i=0; i<2; i++) {
					for(int j=0; j<2; j++) {
						for(int k=0; k<2; k++) {
							// TODO: add an actual interpolation
							res[i][j][k].rho += p.q;
							res[i][j][k].J += p.q;
						}
					}
				}
			}

			// write contributions to contributions grid
			for(int i=0; i<2; i++) {
				for(int j=0; j<2; j++) {
					for(int k=0; k<2; k++) {
						Coord cur = (pos * 2) + Coord{i,j,k};
						contributions[cur] = res[i][j][k];
					}
				}
			}
		}

		/**
		 * Interpolation of particles to grid
		 *
		 * Requires this cell, being located at the position `pos`, to project the effect
		 * of its contained particles to the density grid. Contributions are stored
		 * within the given contributions grid
		 */
		void interP2G(const Coord& pos, allscale::api::user::data::Grid<DensityCell,3>& contributions) const {

			// quick-check
			if (particles.empty()) return;		// nothing to contribute

			// init aggregated densities of neighboring cells
			DensityCell res[2][2][2];
			for(int i=0; i<2; i++) {
				for(int j=0; j<2; j++) {
					for(int k=0; k<2; k++) {
						res[i][j][k].rho = 0.0;
						res[i][j][k].J = 0.0;
					}
				}
			}

			// aggregate particles
			allscale::api::user::pfor(particles, [&](const Particle& p) {
				for(int i=0; i<2; i++) {
					for(int j=0; j<2; j++) {
						for(int k=0; k<2; k++) {
							// TODO: add an actual interpolation
							res[i][j][k].rho += p.q;
							res[i][j][k].J += p.q;
						}
					}
				}
			});

			// write contributions to contributions grid
			for(int i=0; i<2; i++) {
				for(int j=0; j<2; j++) {
					for(int k=0; k<2; k++) {
						Coord cur = (pos * 2)  + Coord{i,j,k};
						contributions[cur] = res[i][j][k];
					}
				}
			}
		}

		/**
		 * This method is updating the position of all particles within this cell for a single
		 * time step, thereby considering the given field as a driving force. Particles
		 * leaving the cell are submitted via channels to neighboring cells.
		 *
		 * @param pos the coordinates of this cell in the grid
		 * @param field the most recently computed state of the surrounding force fields
		 * @param transfers a grid of buffers to send particles to
		 * @param dt the time step to move the particle forward for
		 */
		void moveParticles(const Coord& pos, const Field& field, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers, double dt) {

			// quick-check
			if (particles.empty()) return;

			// -- move the particles in space --

			// extract forces
			Vector<double> E[2][2][2];
			for(int i=0; i<2; i++) {
				for(int j=0; j<2; j++) {
					for(int k=0; k<2; k++) {
						Coord cur({pos[0]+i,pos[1]+j,pos[2]+k});
						E[i][j][k] = field[cur].E;
					}
				}
			}

			// update particles
			allscale::api::user::pfor(particles, [&](Particle& p){

				// TODO: move the computation of forces in an extra function (for unit testing)

				// compute forces
				Vector3<double> f;
				f.x = f.y = f.z = 0.0;
				for(int i=0; i<2; i++) {
					for(int j=0; j<2; j++) {
						for(int k=0; k<2; k++) {
							f.x += E[i][j][k].x * p.q;
							f.y += E[i][j][k].y * p.q;
							f.z += E[i][j][k].z * p.q;
						}
					}
				}

				// update position
				p.updatePosition(dt);

				// update speed
				p.updateVelocity(f,dt);


			});


			// -- migrate particles to other cells if boundaries are crossed --

			// get buffers for particles to be send to neighbors
			Coord size = transfers.size();
			Coord center = pos * 3 + Coord{1,1,1};
			utils::grid<std::vector<Particle>*,3,3,3> neighbors;
			for(int i = 0; i<3; i++) {
				for(int j = 0; j<3; j++) {
					for(int k = 0; k<3; k++) {
						auto cur = center + Coord{i-1,j-1,k-1} * 2;
						// TODO: deal with boundaries
						if (cur[0] < 0 || cur[0] >= size[0]) { neighbors[{i,j,k}] = nullptr; continue; }
						if (cur[1] < 0 || cur[1] >= size[1]) { neighbors[{i,j,k}] = nullptr; continue; }
						if (cur[2] < 0 || cur[2] >= size[2]) { neighbors[{i,j,k}] = nullptr; continue; }
						neighbors[{i,j,k}] = &transfers[cur];
					}
				}
			}

			// move particles
			std::vector<Particle> remaining;
			remaining.reserve(particles.size());
			for(const auto& p : particles) {
				// compute relative position
				double rx = p.x - x;
				double ry = p.y - y;
				double rz = p.z - z;
				if ((fabs(rx) > dx/2) || (fabs(ry) > dy/2) || (fabs(rz) > dz/2)) {
					// compute corresponding neighbor cell
					int i = (rx < -dx/2) ? 0 : ( (rx > dx/2) ? 2 : 1 );
					int j = (ry < -dy/2) ? 0 : ( (ry > dy/2) ? 2 : 1 );
					int k = (rz < -dz/2) ? 0 : ( (rz > dz/2) ? 2 : 1 );
					// send to neighbor cell
					auto target = neighbors[{i,j,k}];
					if (target) target->push_back(p);
				} else {
					// keep particle
					remaining.push_back(p);
				}
			}

			// update content
			particles.swap(remaining);

		}

		/**
 		 * Initial version of the Field Solver: compute fields E and B for the Boris mover
 		 *
 		 * Fields are computed with respect to the center of each cell
 		 */
		void computeFields(const Coord& pos, Vector<double> &E, Vector<double> &B, const bool isDipole){
			if (isDipole) {
				E.x = 0.0;
				E.y = 0.0;
				E.z = 0.0;
				double fac1 = -B0 * pow(Re, 3.0) / pow(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2], 2.5);
				B.x = 3.0 * pos[0] * pos[2] * fac1;
				B.y = 3.0 * pos[1] * pos[2] * fac1;
				B.z = (2.0 * pos[2] * pos[2] - pos[0] * pos[0] - pos[1] * pos[1]) * fac1;
			} else {
				E.x = sin(2.0 * M_PI * pos[0]) * cos(2.0 * M_PI * pos[1]);
				E.y = pos[0] * (1.0 - pos[0]) * pos[1] * (1.0 - pos[1]);
				E.z = pos[0] * pos[0] + pos[2] * pos[2];
				B.x = 0.0;
				B.y = cos(2.0 * M_PI * pos[2]);
				B.z = sin(2.0 * M_PI * pos[0]);
			}
		}

		/**
 		 * Initial version of the Field Solver: compute fields E and B for the Boris mover
 		 *
 		 * Fields are computed with respect to each particle position
 		 */
		void computeFields(const Particle& p, Vector<double> &E, Vector<double> &B, const bool isDipole){
			if (isDipole) {
				E.x = 0.0;
				E.y = 0.0;
				E.z = 0.0;
				double fac1 = -B0 * pow(Re, 3.0) / pow(p.x * p.x + p.y * p.y + p.z * p.z, 2.5);
				B.x = 3.0 * p.x * p.z * fac1;
				B.y = 3.0 * p.y * p.z * fac1;
				B.z = (2.0 * p.z * p.z - p.x * p.x - p.y * p.y) * fac1;
			} else {
				E.x = sin(2.0 * M_PI * p.x) * cos(2.0 * M_PI * p.y);
 				E.y = p.x * (1.0 - p.x) * p.y * (1.0 - p.y);
 				E.z = p.x * p.x + p.z * p.z;
 				B.x = 0.0;
	 			B.y = cos(2.0 * M_PI * p.z);
 				B.z = sin(2.0 * M_PI * p.x);
			}
		}

		/**
 		 * Initialize particles for
 		 * 	isDipole == true: the Earth's dipole simulation
 		 * 	Otherwise:	  the particles and waves interactions
 		 * This step is required for the Boris mover
 		 *
		 * @param time step
		 */
		void initParticles(const Coord&, const double dt, const bool isDipole) {

			// quick-check
			if (particles.empty())
				return;

			if (isDipole) {
				// update particles
				allscale::api::user::pfor(particles, [&](Particle& p){
					Vector<double> v, vr;
					double B_sq, f1, f2;
					double qdto4mc = p.q * dt * 0.25;

					p.q = e; // positive charge, change to -e when simulating electron
					p.q = p.q / m;
					// Trajectory of a proton with 10MeV kinetic energy in dipole field
					//K = K * e;   // convert to Joule
					// Find corresponding speed
					double v_mod = c / sqrt(1.0 + (m * c * c) / K);

					// initial position: equatorial plane 4Re from Earth
					p.x += 4 * Re; p.y += 0.0; p.z += 0.0;
					//p.x = xx[p.id]; p.y = yy[p.id]; p.z = zz[p.id];

					double pitch_angle = 30.0; // initial angle between velocity and mag.field (degrees)
					p.dx = 0.0;
					p.dy = v_mod * sin(pitch_angle * M_PI / 180.0);
					p.dz = v_mod * cos(pitch_angle * M_PI / 180.0);

				    // compute forces
					Vector<double> E, B;
					computeFields(p, E, B, true);
					//computeFields(pos, E, B, true);

					B_sq = B.x * B.x + B.y * B.y + B.z * B.z;
					f1 = tan(qdto4mc * sqrt(B_sq)) / sqrt(B_sq);
					f2 = 2.0 * f1 / (1.0 + f1 * f1 * B_sq);

					// update velocity
					v.x = p.dx + E.x * qdto4mc;
					v.y = p.dy + E.y * qdto4mc;
					v.z = p.dz + E.z * qdto4mc;

					vr.x = v.x + f1 * (v.y * B.z - B.y * v.z);
					vr.y = v.y + f1 * (-v.x * B.z + v.z * B.x);
					vr.z = v.z + f1 * (v.x * B.y - v.y * B.x);
					v.x = v.x + f2 * (vr.y * B.z - vr.z * B.y);
					v.y = v.y + f2 * (-vr.x * B.z + vr.z * B.x);
					v.z = v.z + f2 * (vr.z * B.y - vr.y * B.x);

					p.vxstar = v.x + E.x * qdto4mc;
					p.vystar = v.y + E.y * qdto4mc;
					p.vzstar = v.z + E.z * qdto4mc;
				});
			} else {
				allscale::api::user::pfor(particles, [&](Particle& p){
					Vector<double> v, vr;
					double B_sq, f1, f2;
					double qdto4mc = p.q * dt * 0.25;

					Vector<double> E, B;
					computeFields(p, E, B, false);

					B_sq = B.x * B.x + B.y * B.y + B.z * B.z;
					f1 = tan(qdto4mc * sqrt(B_sq)) / sqrt(B_sq);
					f2 = 2.0 * f1 / (1.0 + f1 * f1 * B_sq);

					// update velocity
					v.x = p.dx + E.x * qdto4mc;
					v.y = p.dy + E.y * qdto4mc;
					v.z = p.dz + E.z * qdto4mc;

					vr.x = v.x + f1 * (v.y * B.z - B.y * v.z);
					vr.y = v.y + f1 * (-v.x * B.z + v.z * B.x);
					vr.z = v.z + f1 * (v.x * B.y - v.y * B.x);
					v.x = v.x + f2 * (vr.y * B.z - vr.z * B.y);
					v.y = v.y + f2 * (-vr.x * B.z + vr.z * B.x);
					v.z = v.z + f2 * (vr.z * B.y - vr.y * B.x);

					p.vxstar = v.x + E.x * qdto4mc;
					p.vystar = v.y + E.y * qdto4mc;
					p.vzstar = v.z + E.z * qdto4mc;
				});
			}
		}

		/**
 		 * Boris mover in cartesian grid
 		 *
		 * This method is updating the position of all particles within this cell for a single
		 * time step, thereby considering the given field as a driving force. Particles
		 * leaving the cell are submitted via channels to neighboring cells.
		 *
		 * @param pos the coordinates of this cell in the grid
		 * @param field the most recently computed state of the surrounding force fields
		 * @param transfers a grid of buffers providing connections to other cells
		 * @param dt a time step
		 * @param isDipole to check with test case we deal with
		 */
		void BorisMover(const Coord& pos, const Field&, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers, const double dt, const bool isDipole) {

			// quick-check
			if (particles.empty())
				return;

			// update particles
			allscale::api::user::pfor(particles, [&](Particle& p){
				Vector<double> v, vr;
				double B_sq, f1, f2;
				double qdto2mc = p.q * dt * 0.5;

				// move particle
				p.x += p.vxstar * dt;
				p.y += p.vystar * dt;
				p.z += p.vzstar * dt;

				// TODO: should not that be a field solver?
				Vector<double> E, B;
				computeFields(p, E, B, isDipole);
				//computeFields(pos, E, B, isDipole);

				B_sq = B.x * B.x + B.y * B.y + B.z * B.z;
				f1 = tan(qdto2mc * sqrt(B_sq)) / sqrt(B_sq);
				f2 = 2.0 * f1 / (1.0 + f1 * f1 * B_sq);

				// update velocity
				v.x = p.vxstar + E.x * qdto2mc;
				v.y = p.vystar + E.y * qdto2mc;
				v.z = p.vzstar + E.z * qdto2mc;

				vr.x = v.x + f1 * (v.y * B.z - B.y * v.z);
				vr.y = v.y + f1 * (-v.x * B.z + v.z * B.x);
				vr.z = v.z + f1 * (v.x * B.y - v.y * B.x);
				v.x = v.x + f2 * (vr.y * B.z - vr.z * B.y);
				v.y = v.y + f2 * (-vr.x * B.z + vr.z * B.x);
				v.z = v.z + f2 * (vr.x * B.y - vr.y * B.x);

				vr.x = v.x + E.x * qdto2mc;
				vr.y = v.y + E.y * qdto2mc;
				vr.z = v.z + E.z * qdto2mc;

				p.dx = (p.vxstar + vr.x) * 0.5;
				p.dy = (p.vystar + vr.y) * 0.5;
				p.dz = (p.vzstar + vr.z) * 0.5;

				p.vxstar = vr.x;
				p.vystar = vr.y;
				p.vzstar = vr.z;
			});

			// TODO: test this particles exchanger
			// 		161208: seems to work
			// get buffers for particles to be send to neighbors
			Coord size = transfers.size();
			utils::grid<std::vector<Particle>*,3,3,3> neighbors;
			Coord center = pos * 3 + Coord{1,1,1};
			for(int i = 0; i<3; i++) {
				for(int j = 0; j<3; j++) {
					for(int k = 0; k<3; k++) {
						neighbors[{i,j,k}] = nullptr;
						auto cur = center + Coord{i-1,j-1,k-1} * 2;
						if (cur[0] < 0 || cur[0] >= size[0]) continue;
						if (cur[1] < 0 || cur[1] >= size[1]) continue;
						if (cur[2] < 0 || cur[2] >= size[2]) continue;
						neighbors[{i,j,k}] = &transfers[cur];
					}
				}
			}

			// move particles
			std::vector<Particle> remaining;
			remaining.reserve(particles.size());
			for(const auto& p : particles) {
				// compute relative position
				double rx = p.x - x;
				double ry = p.y - y;
				double rz = p.z - z;

				if ((fabs(rx) > dx/2) || (fabs(ry) > dy/2) || (fabs(rz) > dz/2)) {
					// compute corresponding neighbor cell
					int i = (rx < 0) ? 0 : 2;
					int j = (ry < 0) ? 0 : 2;
					int k = (rz < 0) ? 0 : 2;

					// send to neighbor cell
					auto target = neighbors[{i,j,k}];
					if (target) target->push_back(p);

				} else {
					// keep particle
					remaining.push_back(p);
				}
			}

			// update content
			// TODO: what does this swap do?
			particles.swap(remaining);
		}

		/**
 		 * Interpolation of fields to particles and the Boris mover in cartesian grid
 		 *
		 * This method is updating the position of all particles within this cell for a single
		 * time step, thereby considering the given field as a driving force. Particles
		 * leaving the cell are submitted via channels to neighboring cells.
		 *
		 * @param pos the coordinates of this cell in the grid
		 * @param field the most recently computed state of the surrounding force fields
		 * @param transfers a grid of buffers to send particles to
		 * @param time step
		 */
		void InterF2PBorisMover(const Coord& pos, const Field& field, allscale::api::user::data::Grid<std::vector<Particle>,3>&, const double dt, const bool isDipole) {

			// quick-check
			if (particles.empty())
				return;

			// extract forces
			Vector<double> E[2][2][2];
			Vector<double> B[2][2][2];
			for(int i=0; i<2; i++) {
				for(int j=0; j<2; j++) {
					for(int k=0; k<2; k++) {
						Coord cur({pos[0]+i,pos[1]+j,pos[2]+k});
						E[i][j][k] = field[cur].E;
						B[i][j][k] = field[cur].B;
					}
				}
			}

			// update particles
			allscale::api::user::pfor(particles, [&](Particle& p){
				Vector<double> v, vr;
				double B_sq, f1, f2;
				const double qdto2mc = p.q * dt * 0.5;

				// BEGIN of interpolation of fields to particles
				Vector<double> corig, vorig, cavg, vavg;
				corig.x = p.x;
				corig.y = p.y;
				corig.z = p.z;
				vorig.x = p.dx;
				vorig.y = p.dy;
				vorig.z = p.dz;
				cavg = corig;
				vavg = vorig;

				// END of interpolation of fields to particles

				// move particle
				p.x += p.vxstar * dt;
				p.y += p.vystar * dt;
				p.z += p.vzstar * dt;

				Vector<double> E, B;
				computeFields(p, E, B, isDipole);

				B_sq = B.x * B.x + B.y * B.y + B.z * B.z;
				f1 = tan(qdto2mc * sqrt(B_sq)) / sqrt(B_sq);
				f2 = 2.0 * f1 / (1.0 + f1 * f1 * B_sq);

				// update velocity
				v.x = p.vxstar + E.x * qdto2mc;
				v.y = p.vystar + E.y * qdto2mc;
				v.z = p.vzstar + E.z * qdto2mc;

				vr.x = v.x + f1 * (v.y * B.z - B.y * v.z);
				vr.y = v.y + f1 * (-v.x * B.z + v.z * B.x);
				vr.z = v.z + f1 * (v.x * B.y - v.y * B.x);
				v.x = v.x + f2 * (vr.y * B.z - vr.z * B.y);
				v.y = v.y + f2 * (-vr.x * B.z + vr.z * B.x);
				v.z = v.z + f2 * (vr.x * B.y - vr.y * B.x);

				vr.x = v.x + E.x * qdto2mc;
				vr.y = v.y + E.y * qdto2mc;
				vr.z = v.z + E.z * qdto2mc;

				p.dx = (p.vxstar + vr.x) * 0.5;
				p.dy = (p.vystar + vr.y) * 0.5;
				p.dz = (p.vzstar + vr.z) * 0.5;

				p.vxstar = vr.x;
				p.vystar = vr.y;
				p.vzstar = vr.z;
			});

			// TODO: test this particles exchanger
			/*// get buffers for particles to be send to neighbors
			Coord size = transfers.getSize();
			utils::grid<std::vector<Particle>*,3,3,3> neighbors;
			Coord center = pos * 3 + Coord{1,1,1};
			for(int i = 0; i<3; i++) {
				for(int j = 0; j<3; j++) {
					for(int k = 0; k<3; k++) {
						auto cur = center + Coord{i-1,j-1,k-1} * 2;
						if (cur[0] < 0 || cur[0] >= size[0]) continue;
						if (cur[1] < 0 || cur[1] >= size[1]) continue;
						if (cur[2] < 0 || cur[2] >= size[2]) continue;
						neighbors[{i,j,k}] = &transfers[cur];
					}
				}
			}

			// move particles
			std::vector<Particle> remaining;
			remaining.reserve(particles.size());
			for(const auto& p : particles) {

				// compute relative position
				double rx = p.x - x;
				double ry = p.y - y;
				double rz = p.z - z;

				if ((fabs(rx) > dx/2) || (fabs(ry) > dy/2) || (fabs(rz) > dz/2)) {

					// compute corresponding neighbor cell
					int i = (rx < 0) ? 0 : 2;
					int j = (ry < 0) ? 0 : 2;
					int k = (rz < 0) ? 0 : 2;

					// send to neighbor cell
					auto target = neighbors[{i,j,k}];
					if (target) target->push_back(p);

				} else {
					// keep particle
					remaining.push_back(p);
				}
			}

			// update content
			particles.swap(remaining);*/
		}

		/**
		 * Imports the particles from directed towards this cell from the given transfere buffer.
		 */
		void importParticles(const Coord& pos, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers) {

			// import particles send to this cell
			Coord size = transfers.size();
			Coord center = pos * 3 + Coord{1,1,1};
			for(int i = 0; i<3; i++) {
				for(int j = 0; j<3; j++) {
					for(int k = 0; k<3; k++) {
						auto cur = center + Coord{i-1,j-1,k-1};
						if (cur[0] < 0 || cur[0] >= size[0]) continue;
						if (cur[1] < 0 || cur[1] >= size[1]) continue;
						if (cur[2] < 0 || cur[2] >= size[2]) continue;
						auto& in = transfers[cur];
						particles.insert(particles.end(), in.begin(), in.end());
						in.clear();
					}
				}
			}

		}

	};

	using Cells = allscale::api::user::data::Grid<Cell,3>;



} // end namespace ipic3d
