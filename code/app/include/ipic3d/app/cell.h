#pragma once

#include <vector>

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/data/vector.h"
#include "allscale/api/user/operator/pfor.h"

#include "ipic3d/app/particle.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/universe_properties.h"
#include "ipic3d/app/utils/points.h"
#include "ipic3d/app/utils/static_grid.h"
#include "ipic3d/app/parameters.h"

namespace ipic3d {


	struct DensityCell {
		double rho;			// charge density
		Vector3<double> J;	// current density
	};

	using Density = allscale::api::user::data::Grid<DensityCell,3>;	// a 3D grid of density cells

	/**
	 * The structure of a single cell, forming a container for a set of particles
	 * located within a confined area of space.
	 */
	struct Cell {

		using Coord = utils::Coordinate<3>;

		// the list of local particles
		std::vector<Particle> particles;

		/**
		 * Requires this cell, being located at position `pos`, to project the effect
		 * of its contained particles to the density grid. Contributions are stored
		 * within the given contributions grid
		 */
		void projectToDensityField(const UniverseProperties& universeProperties, const Coord& pos, DensityCell& contributions) const {

			// quick-check
			if (particles.empty()) return;		// nothing to contribute

			// init aggregated densities of neighboring cells
			contributions.rho = 0.0;
			contributions.J = 0.0;

			// aggregate particles
			// TODO data race on res, should be avoided
			for(const auto& p : particles) {
				contributions.rho += p.q;

				contributions.J.x += p.q * p.velocity.x;
				contributions.J.y += p.q * p.velocity.y;
				contributions.J.z += p.q * p.velocity.z;
			}

        	double vol = universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;
			contributions.rho = contributions.rho / vol;
			contributions.J.x = contributions.J.x / vol;
			contributions.J.y = contributions.J.y / vol;
			contributions.J.z = contributions.J.z / vol;
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
		void initParticles(const UniverseProperties& /*universeProperties*/, const Coord&, const double dt, const bool isDipole) {

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
		void BorisMover(const UniverseProperties& universeProperties, const Coord& pos, const Field&, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers, const double dt, const bool isDipole) {

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
				Vector3<double> relPos = p.position - getCenterOfCell(pos, universeProperties);
				auto halfWidth = universeProperties.cellWidth/2;
				if ((fabs(relPos.x) > halfWidth.x) || (fabs(relPos.y) > halfWidth.y) || (fabs(relPos.z) > halfWidth.z)) {
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
	void moveParticlesFirstOrder(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field) {

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
			p.updatePosition(universeProperties.dt);

			// update speed
			p.updateVelocity(f,universeProperties.dt);


		});

	}

	/**
	* This method computes a trilinear interpolation for a given target position within a rectangular box spanned by 8 corner points
	*   This is interpolation of fields to particles
	* Math: http://paulbourke.net/miscellaneous/interpolation/
	* TODO: Move this to some math utilities header?
	*
	* @param corners the 8 surrounding points to interpolate from
	* @param pos the target position for which to interpolate
	*/
	template<typename T>
	T trilinearInterpolationF2P(const T corners[2][2][2], const Vector3<double>& pos, const double vol) {
		T res = T();

	    for(int i = 0; i < 2; ++i) {
		    for(int j = 0; j < 2; ++j) {
			    for(int k = 0; k < 2; ++k) {
				    auto fac = (i == 0 ? (1 - pos.x) : pos.x) * (j == 0 ? (1 - pos.y) : pos.y) * (k == 0 ? (1 - pos.z) : pos.z) / vol;
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
	void moveParticlesBorisStyle(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field) {

		// quick-check
		if (cell.particles.empty()) return;

		// -- move the particles in space --

		// extract forces
		// TODO: move this to some C++ structure
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

		const auto cellCenter = getCenterOfCell(pos, universeProperties);

        double vol = 1.0; //universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;
		// update particles
		allscale::api::user::pfor(cell.particles, [&](Particle& p){

			// Docu: https://www.particleincell.com/2011/vxb-rotation/
			// Code: https://www.particleincell.com/wp-content/uploads/2011/07/ParticleIntegrator.java

			// get relative position of particle within cell
			const auto relPos = allscale::api::user::data::elementwiseDivision((p.position - (cellCenter - universeProperties.cellWidth*0.5)), (universeProperties.cellWidth));

			// interpolate
			auto E = trilinearInterpolationF2P(Es, relPos, vol);
			auto B = trilinearInterpolationF2P(Bs, relPos, vol);

			// update velocity
			p.updateVelocityBorisStyle(E, B, universeProperties.dt);

			// update position
			p.updatePosition(universeProperties.dt);

		});

	}


	/**
	 * This method extracts all particles which are no longer in the domain of the
	 * given cell and inserts them into the provided transfer buffers.
	 *
	 * @param pos the coordinates of this cell in the grid
	 * @param transfers a grid of buffers to send particles to
	 */
	void exportParticles(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers) {

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
			Vector3<double> relPos = p.position - getCenterOfCell(pos, universeProperties);
			auto halfWidth = universeProperties.cellWidth / 2;
			if ((fabs(relPos.x) > halfWidth.x) || (fabs(relPos.y) > halfWidth.y) || (fabs(relPos.z) > halfWidth.z)) {
				// compute corresponding neighbor cell
				int i = (relPos.x < -halfWidth.x) ? 0 : ( (relPos.x > halfWidth.x) ? 2 : 1 );
				int j = (relPos.y < -halfWidth.y) ? 0 : ( (relPos.y > halfWidth.y) ? 2 : 1 );
				int k = (relPos.z < -halfWidth.z) ? 0 : ( (relPos.z > halfWidth.z) ? 2 : 1 );
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
	void importParticles(const UniverseProperties& /*universeProperties*/, Cell& cell, const utils::Coordinate<3>& pos, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers) {

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


	/**
	 * Static Field Solver: Fields are computed with respect to the center of each cell
	 */
	void FieldSolverStatic(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, Field& field) {
		switch(universeProperties.useCase) {

			case UseCase::Dipole:
				{
					field[pos].E = {0.0, 0.0, 0.0};
					auto B = field[pos].B;
					auto cellCenter = getCenterOfCell(pos, universeProperties);

					double fac1 = -B0 * pow(Re, 3.0) / pow(allscale::api::user::data::sumOfSquares(cellCenter), 2.5);
					B.x = 3.0 * cellCenter.x * cellCenter.z * fac1;
					B.y = 3.0 * cellCenter.y * cellCenter.z * fac1;
					B.z = (2.0 * cellCenter.z * cellCenter.z - cellCenter.x * cellCenter.x - cellCenter.y * cellCenter.y) * fac1;

					field[pos].B = B;
				}
				 break;

			case UseCase::ParticleWave:
				{
					// TODO: to provide
					break;
				}

			case UseCase::Test:
				{
					// TODO: to provide
					break;
				}

			default:
				assert_not_implemented() << "The specified use case is not supported yet!";
		}
	}

    /**
     * Explicit Field Solver: Fields are computed using leapfrog algorithm
     */
	void FieldSolverLeapfrog(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, Field& field) {
        // 1. Compute current density J as sum of particles density times particles velocity
        
        // 2. Compute electic field E using leapfrog with the time step delta t 
        
        // 3. Compute magnetic field B using leapfrog with the time step delta t, but starts on detla t / 2
        //    Compute also magnetic field B on the center of each cell as average of all nodes
    }

} // end namespace ipic3d
