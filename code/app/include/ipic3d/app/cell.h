#pragma once

#include <vector>

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/data/vector.h"
#include "allscale/api/user/operator/pfor.h"

#include "ipic3d/app/particle.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/utils/points.h"
#include "ipic3d/app/utils/static_grid.h"

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

		// center of the cell
		Vector3<double> center;

		// the cell grid spacing
		Vector3<double> spacing;

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
 		 * Initial version of the Field Solver: compute fields E and B for the Boris mover
 		 *
 		 * Fields are computed with respect to the center of each cell
 		 */
	    void computeFields(const Coord& pos, Vector3<double>& E, Vector3<double>& B, const bool isDipole) {
		    if(isDipole) {
			    E = {0, 0, 0};
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
		void computeFields(const Particle& p, Vector3<double> &E, Vector3<double> &B, const bool isDipole){
			if (isDipole) {
			    E = {0, 0, 0};
			    double fac1 = -B0 * pow(Re, 3.0) / pow(sumOfSquares(p.position), 2.5);
				B.x = 3.0 * p.position.x * p.position.z * fac1;
				B.y = 3.0 * p.position.y * p.position.z * fac1;
				B.z = (2.0 * p.position.z * p.position.z - p.position.x * p.position.x - p.position.y * p.position.y) * fac1;
			} else {
				E.x = sin(2.0 * M_PI * p.position.x) * cos(2.0 * M_PI * p.position.y);
 				E.y = p.position.x * (1.0 - p.position.x) * p.position.y * (1.0 - p.position.y);
 				E.z = p.position.x * p.position.x + p.position.z * p.position.z;
 				B.x = 0.0;
	 			B.y = cos(2.0 * M_PI * p.position.z);
 				B.z = sin(2.0 * M_PI * p.position.x);
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
					Vector3<double> v, vr;
					double B_sq, f1, f2;
					double qdto4mc = p.q * dt * 0.25;

					p.q = e; // positive charge, change to -e when simulating electron
					p.q = p.q / m;
					// Trajectory of a proton with 10MeV kinetic energy in dipole field
					//K = K * e;   // convert to Joule
					// Find corresponding speed
					double v_mod = c / sqrt(1.0 + (m * c * c) / K);

					// initial position: equatorial plane 4Re from Earth
					p.position.x += 4 * Re; p.position.y += 0.0; p.position.z += 0.0;

					double pitch_angle = 30.0; // initial angle between velocity and mag.field (degrees)
					p.velocity.x = 0.0;
					p.velocity.y = v_mod * sin(pitch_angle * M_PI / 180.0);
					p.velocity.z = v_mod * cos(pitch_angle * M_PI / 180.0);

				    // compute forces
					Vector3<double> E, B;
					computeFields(p, E, B, true);
					//computeFields(pos, E, B, true);

					B_sq = sumOfSquares(B);
					f1 = tan(qdto4mc * sqrt(B_sq)) / sqrt(B_sq);
					f2 = 2.0 * f1 / (1.0 + f1 * f1 * B_sq);

					// update velocity
					v = p.velocity + E * qdto4mc;
					vr = v + f1 * crossProduct(v, B);
					v = v + f2 * crossProduct(vr, B);

					p.velocityStar = v + E * qdto4mc;
				});
			} else {
				allscale::api::user::pfor(particles, [&](Particle& p){
					Vector3<double> v, vr;
					double B_sq, f1, f2;
					double qdto4mc = p.q * dt * 0.25;

					Vector3<double> E, B;
					computeFields(p, E, B, false);

					B_sq = sumOfSquares(B);
					f1 = tan(qdto4mc * sqrt(B_sq)) / sqrt(B_sq);
					f2 = 2.0 * f1 / (1.0 + f1 * f1 * B_sq);

					// update velocity
					v = p.velocity + E * qdto4mc;

					vr = v + f1 * crossProduct(v, B);
					v = v + f2 * crossProduct(vr, B);

					p.velocityStar = v + E * qdto4mc;
				});
			}
		}

		/**
 		 * Boris mover in Cartesian grid
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
				Vector3<double> v, vr;
				double B_sq, f1, f2;
				double qdto2mc = p.q * dt * 0.5;

				// move particle
				p.position += p.velocityStar * dt;

				// TODO: should not that be a field solver?
				Vector3<double> E, B;
				computeFields(p, E, B, isDipole);
				//computeFields(pos, E, B, isDipole);

				B_sq = sumOfSquares(B);
				f1 = tan(qdto2mc * sqrt(B_sq)) / sqrt(B_sq);
				f2 = 2.0 * f1 / (1.0 + f1 * f1 * B_sq);

				// update velocity
				v = p.velocityStar + E * qdto2mc;

				vr = v + f1 * crossProduct(v, B);
				v = v + f2 * crossProduct(vr, B);

			    vr = v + E * qdto2mc;

				p.velocity = (p.velocityStar + vr) * 0.5;

				p.velocityStar = vr;
			});

			// get buffers for particles to be send to neighbors
			Coord size = transfers.size();
			utils::grid<std::vector<Particle>*,3,3,3> neighbors;
			Coord centerIndex = pos * 3 + Coord{1,1,1};
			for(int i = 0; i<3; i++) {
				for(int j = 0; j<3; j++) {
					for(int k = 0; k<3; k++) {
						neighbors[{i,j,k}] = nullptr;
						auto cur = centerIndex + Coord{i-1,j-1,k-1} * 2;
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
				Vector3<double> relPos = p.position - center;

				if ((fabs(relPos.x) > spacing.x/2) || (fabs(relPos.y) > spacing.y/2) || (fabs(relPos.z) > spacing.z/2)) {
					// compute corresponding neighbor cell
					int i = (relPos.x < 0) ? 0 : 2;
					int j = (relPos.y < 0) ? 0 : 2;
					int k = (relPos.z < 0) ? 0 : 2;

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
 		 * Interpolation of fields to particles and the Boris mover in Cartesian grid
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
			Vector3<double> E[2][2][2];
			Vector3<double> B[2][2][2];
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
				Vector3<double> v, vr;
				double B_sq, f1, f2;
				const double qdto2mc = p.q * dt * 0.5;

				// move particle
				p.position += p.velocityStar * dt;

				Vector3<double> E, B;
				computeFields(p, E, B, isDipole);

				B_sq = sumOfSquares(B);
				f1 = tan(qdto2mc * sqrt(B_sq)) / sqrt(B_sq);
				f2 = 2.0 * f1 / (1.0 + f1 * f1 * B_sq);

				// update velocity
				v = p.velocityStar + E * qdto2mc;

				vr = v + f1 * crossProduct(v, B);
				v = v + f2 * crossProduct(vr, B);

				vr = v + E * qdto2mc;

				p.velocity = (p.velocityStar + vr) * 0.5;

				p.velocityStar = vr;
			});

		}

	};

	using Cells = allscale::api::user::data::Grid<Cell,3>;


	/**
	 * This method is updating the position of all particles within a cell for a single
	 * time step, thereby considering the given field as a driving force.
	 *
	 * @param pos the coordinates of this cell in the grid
	 * @param field the most recently computed state of the surrounding force fields
	 * @param dt the time step to move the particle forward for
	 */
	void moveParticlesFirstOrder(Cell& cell, const utils::Coordinate<3>& pos, const Field& field, double dt) {

		// quick-check
		if (cell.particles.empty()) return;

		// -- move the particles in space --

		// extract forces
		Vector3<double> E[2][2][2];
		for(int i=0; i<2; i++) {
			for(int j=0; j<2; j++) {
				for(int k=0; k<2; k++) {
					utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});
					E[i][j][k] = field[cur].E;
				}
			}
		}

		// update particles
		allscale::api::user::pfor(cell.particles, [&](Particle& p){

			// TODO: move the computation of forces in an extra function (for unit testing)

			// compute forces
			Vector3<double> f;
			f.x = f.y = f.z = 0.0;
			for(int i=0; i<2; i++) {
				for(int j=0; j<2; j++) {
					for(int k=0; k<2; k++) {
						f += E[i][j][k] * p.q;
					}
				}
			}

			// update position
			p.updatePosition(dt);

			// update speed
			p.updateVelocity(f,dt);


		});

	}

	/**
	* This method computes a trilinear interpolation for a given target position within a rectangular box spanned by 8 corner points
	* Math: http://paulbourke.net/miscellaneous/interpolation/
	* TODO: Move this to some math utilities header?
	*
	* @param corners the 8 surrounding points to interpolate from
	* @param pos the target position for which to interpolate
	*/
	template<typename T>
	T trilinearInterpolation(const T corners[2][2][2], const Vector3<double>& pos) {
		T res = T();

	    for(int i = 0; i < 2; ++i) {
		    for(int j = 0; j < 2; ++j) {
			    for(int k = 0; k < 2; ++k) {
				    auto fac = (i == 0 ? (1 - pos.x) : pos.x) * (j == 0 ? (1 - pos.y) : pos.y) * (k == 0 ? (1 - pos.z) : pos.z);
				    res += corners[i][j][k] * fac;
			    }
		    }
	    }

	    return res;
	}

	/**
	 * This method is updating the position of all particles within a cell for a single
	 * time step, thereby considering the given field as a driving force.
	 *
	 * @param pos the coordinates of this cell in the grid
	 * @param field the most recently computed state of the surrounding force fields
	 * @param dt the time step to move the particle forward for
	 */
	void moveParticlesBorisStyle(Cell& cell, const utils::Coordinate<3>& pos, const Field& field, double dt) {

		// quick-check
		if (cell.particles.empty()) return;

		// -- move the particles in space --

		// extract forces
		Vector3<double> Es[2][2][2];
		Vector3<double> Bs[2][2][2];
		for(int i=0; i<2; i++) {
			for(int j=0; j<2; j++) {
				for(int k=0; k<2; k++) {
					utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});
					Es[i][j][k] = field[cur].E;
					Bs[i][j][k] = field[cur].B;
				}
			}
		}

		// update particles
		allscale::api::user::pfor(cell.particles, [&](Particle& p){

			// Docu: https://www.particleincell.com/2011/vxb-rotation/
			// Code: https://www.particleincell.com/wp-content/uploads/2011/07/ParticleIntegrator.java

			// get relative position of particle within cell
			const auto relPos = allscale::api::user::data::elementwiseDivision((p.position - (cell.center - cell.spacing*0.5)), (cell.spacing));

			// interpolate
			auto E = trilinearInterpolation(Es, relPos);
			auto B = trilinearInterpolation(Bs, relPos);

			// update velocity
			p.updateVelocityBorisStyle(E, B, dt);

			// update position
			p.updatePosition(dt);

		});

	}


	/**
	 * This method extracts all particles which are no longer in the domain of the
	 * given cell and inserts them into the provided transfer buffers.
	 *
	 * @param pos the coordinates of this cell in the grid
	 * @param transfers a grid of buffers to send particles to
	 */
	void exportParticles(Cell& cell, const utils::Coordinate<3>& pos, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers) {


		// -- migrate particles to other cells if boundaries are crossed --

		// get buffers for particles to be send to neighbors
		utils::Coordinate<3> size = transfers.size();
		utils::Coordinate<3> centerIndex = pos * 3 + utils::Coordinate<3>{1,1,1};
		utils::grid<std::vector<Particle>*,3,3,3> neighbors;
		for(int i = 0; i<3; i++) {
			for(int j = 0; j<3; j++) {
				for(int k = 0; k<3; k++) {
					auto cur = centerIndex + utils::Coordinate<3>{i-1,j-1,k-1} * 2;
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
		remaining.reserve(cell.particles.size());
		for(const auto& p : cell.particles) {
			// compute relative position
			Vector3<double> relPos = p.position - cell.center;
			if ((fabs(relPos.x) > cell.spacing.x/2) || (fabs(relPos.y) > cell.spacing.y/2) || (fabs(relPos.z) > cell.spacing.z/2)) {
				// compute corresponding neighbor cell
				int i = (relPos.x < -cell.spacing.x/2) ? 0 : ( (relPos.x > cell.spacing.x/2) ? 2 : 1 );
				int j = (relPos.y < -cell.spacing.y/2) ? 0 : ( (relPos.y > cell.spacing.y/2) ? 2 : 1 );
				int k = (relPos.z < -cell.spacing.z/2) ? 0 : ( (relPos.z > cell.spacing.z/2) ? 2 : 1 );
				// send to neighbor cell
				auto target = neighbors[{i,j,k}];
				if (target) target->push_back(p);
			} else {
				// keep particle
				remaining.push_back(p);
			}
		}

		// update content
		cell.particles.swap(remaining);

	}


	/**
	 * Imports the particles from directed towards this cell from the given transfer buffer.
	 */
	void importParticles(Cell& cell, const utils::Coordinate<3>& pos, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers) {

		// import particles send to this cell
		utils::Coordinate<3> size = transfers.size();
		utils::Coordinate<3> centerIndex = pos * 3 + utils::Coordinate<3>{1,1,1};
		for(int i = 0; i<3; i++) {
			for(int j = 0; j<3; j++) {
				for(int k = 0; k<3; k++) {
					auto cur = centerIndex + utils::Coordinate<3>{i-1,j-1,k-1};
					if (cur[0] < 0 || cur[0] >= size[0]) continue;
					if (cur[1] < 0 || cur[1] >= size[1]) continue;
					if (cur[2] < 0 || cur[2] >= size[2]) continue;
					auto& in = transfers[cur];
					cell.particles.insert(cell.particles.end(), in.begin(), in.end());
					in.clear();
				}
			}
		}

	}



} // end namespace ipic3d
