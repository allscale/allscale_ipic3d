#pragma once

#include <vector>
#include <random>

#include "allscale/api/user/data/grid.h"
#include "allscale/utils/vector.h"
#include "allscale/api/user/operator/pfor.h"

#include "ipic3d/app/particle.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/universe_properties.h"
#include "ipic3d/app/utils/points.h"
#include "ipic3d/app/utils/static_grid.h"
#include "ipic3d/app/parameters.h"

namespace ipic3d {

	//using Coord = utils::Coordinate<3>;

	/**
	 * The structure of a single cell, forming a container for particles
	 * located within a confined (rectangular) area of space.
	 */
	struct Cell {

		// the list of local particles
		std::vector<Particle> particles;

	};

	using Cells = allscale::api::user::data::Grid<Cell, 3>; // a 3D grid of cells


	Cells initCells(const InitProperties& initProperties, const UniverseProperties& properties) {

		const utils::Coordinate<3> zero = 0;							// a zero constant (coordinate [0,0,0])

		// -- initialize the grid of cells --

		// the 3-D grid of cells
		//Grid<Cell> cells(properties.size);								// the grid of cells containing the particles
		Cells cells(properties.size);

		// -- initialize the state of each individual cell --

		// TODO: return this as a treeture
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {

			Cell& cell = cells[pos];

			// -- add particles --

			// compute number of particles to be added
			unsigned particlesPerCell = initProperties.particlesPerCell[0].x + initProperties.particlesPerCell[0].y + initProperties.particlesPerCell[0].z;

			// add the requested number of parameters
			std::minstd_rand randGenerator((unsigned)(pos[0] * 10000 + pos[1] * 100 + pos[2]));
			const double randMax = std::minstd_rand::max();
			for (unsigned i = 0; i < particlesPerCell; i++) {
				Particle p;

				Vector3<double> randVals = {(double)randGenerator() / randMax, (double)randGenerator() / randMax, (double)randGenerator() / randMax};
				// initialize particle position
				p.position = getCenterOfCell(pos, properties) + allscale::utils::elementwiseProduct(randVals, properties.cellWidth) - properties.cellWidth/2;

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

	/**
	* This function projects the effect of the particles contained in the specified cell
	* to the given density contributions grid.
	*
	* @param universeProperties the properties of this universe
	* @param cell the cell whose particles are considered in the density contributions computation
	* @param pos the coordinates of this cell in the grid
	* @param contributions the density contributions output
	*/
	void projectToDensityField(const UniverseProperties& universeProperties, const Cell& cell, const utils::Coordinate<3>& pos, DensityNodes& density) {

		// quick-check
		if(cell.particles.empty()) return;		// nothing to contribute

		// init aggregated densities of neighboring cells
		Vector3<double> Js[2][2][2] = {0};

		// aggregate charge density from particles
		// TODO: use pfor here
		for(const auto& p : cell.particles) {
			// computation of J also includes weights from the particles as for E
			for(int i = 0; i < 2; ++i) {
				for(int j = 0; j < 2; ++j) {
					for(int k = 0; k < 2; ++k) {
						utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});
						auto cornerPos = getLocationForCells(cur, universeProperties); 
						auto fac = (i == 0 ? (1 - cornerPos.x) : cornerPos.x) * (j == 0 ? (1 - cornerPos.y) : cornerPos.y) * (k == 0 ? (1 - cornerPos.z) : cornerPos.z);
						Js[i][j][k] += p.q * p.velocity * fac;
					}
				}
			}
		}

		double vol = universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;
		for(int i=0; i<2; i++) {
			for(int j=0; j<2; j++) {
				for(int k=0; k<2; k++) {
					utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});
					density[cur].J = Js[i][j][k] / vol;
				}
			}
		}
	}

	/**
	 * This function updates the position of all particles within a cell for a single
	 * time step, considering the given field as a driving force.
	 *
	 * @param universeProperties the properties of this universe
	 * @param the cell whose particles are to be moved
	 * @param pos the coordinates of this cell in the grid
	 * @param field the most recently computed state of the surrounding force fields
	 */
	void moveParticlesFirstOrder(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field) {

		assert_true(pos.dominatedBy(universeProperties.size)) << "Position " << pos << " is outside universe of size " << universeProperties.size;

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
		allscale::api::user::pfor(cell.particles, [&](Particle& p) {

			// update position
			p.updatePosition(universeProperties.dt);

			// compute electric force
			auto f = computeElectricForce(E, p);

			// update speed
			p.updateVelocity(f,universeProperties.dt);


		});

	}

	/**
	* This function computes a trilinear interpolation for a given target position within a rectangular box spanned by 8 corner points
	*   This is interpolation of fields to particles
	* Math: http://paulbourke.net/miscellaneous/interpolation/
	* TODO: Move this to some math utilities header?
	*
	* @param corners the 8 surrounding points to interpolate from
	* @param pos the target position for which to interpolate
	* @param vol the volume spanned by the 8 surrounding points
	*/
	template<typename T>
	T trilinearInterpolationF2P(const T corners[2][2][2], const Vector3<double>& pos, const double vol) {
		T res = T(0);

	    for(int i = 0; i < 2; ++i) {
		    for(int j = 0; j < 2; ++j) {
			    for(int k = 0; k < 2; ++k) {
				    auto fac = (i == 0 ? (1 - pos.x) : pos.x) * (j == 0 ? (1 - pos.y) : pos.y) * (k == 0 ? (1 - pos.z) : pos.z);
				    res += corners[i][j][k] * fac;
			    }
		    }
	    }

	    return res / vol;
	}

	/**
	 * This function updates the position of all particles within a cell for a single
	 * time step, considering the given field as a driving force.
	 *
	 * @param universeProperties the properties of this universe
	 * @param cell the cell whose particles are moved
	 * @param pos the coordinates of this cell in the grid
	 * @param field the most recently computed state of the surrounding force fields
	 */
	void moveParticlesBorisStyle(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field) {

		assert_true(pos.dominatedBy(universeProperties.size)) << "Position " << pos << " is outside universe of size " << universeProperties.size;

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
					cur += utils::Coordinate<3>(1); // shift because of the boundary cells
					Es[i][j][k] = field[cur].E;
					Bs[i][j][k] = field[cur].B;
				}
			}
		}

		const auto cellCenter = getCenterOfCell(pos, universeProperties);

        double vol = universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;

		// update particles
		allscale::api::user::pfor(cell.particles, [&](Particle& p){
			// Docu: https://www.particleincell.com/2011/vxb-rotation/
			// Code: https://www.particleincell.com/wp-content/uploads/2011/07/ParticleIntegrator.java

			// get relative position of particle within cell
			const auto relPos = allscale::utils::elementwiseDivision((p.position - (cellCenter - universeProperties.cellWidth*0.5)), (universeProperties.cellWidth));

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
	 * This function extracts all particles which are no longer in the domain of the
	 * given cell and inserts them into the provided transfer buffers.
	 *
	 * @param universeProperties the properties of this universe
	 * @param cell the cell whose particles are moved
	 * @param pos the coordinates of this cell in the grid
	 * @param transfers a grid of buffers to send particles to
	 */
	void exportParticles(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers) {

		assert_true(pos.dominatedBy(universeProperties.size)) << "Position " << pos << " is outside universe of size " << universeProperties.size;

		// -- migrate particles to other cells if boundaries are crossed --

		// get buffers for particles to be send to neighbors
		utils::Coordinate<3> size = transfers.size();
		utils::Coordinate<3> centerIndex = pos * 3 + utils::Coordinate<3>{1,1,1};
		utils::grid<std::vector<Particle>*,3,3,3> neighbors;
		for(int i = 0; i<3; i++) {
			for(int j = 0; j<3; j++) {
				for(int k = 0; k<3; k++) {
					auto cur = centerIndex + utils::Coordinate<3>{i-1,j-1,k-1} * 2;
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
		for(auto& p : cell.particles) {
			// compute relative position
			Vector3<double> relPos = p.position - getCenterOfCell(pos, universeProperties);
			auto halfWidth = universeProperties.cellWidth / 2;
			if ((fabs(relPos.x) > halfWidth.x) || (fabs(relPos.y) > halfWidth.y) || (fabs(relPos.z) > halfWidth.z)) {
				// compute corresponding neighbor cell
				// cover the inner cells as well as the boundary cells on positions 0 and N-1
				int i = (relPos.x < -halfWidth.x) ? ( pos.x == 0 ? size.x - 1 : 0 ) : ( (relPos.x > halfWidth.x) ? ( (pos.x == universeProperties.size.x - 1) ? 0 : 2 ) : 1 );
				// adjust particle's position in case it exits the domain
				p.position.x = ( (pos.x == 0) && (relPos.x < -halfWidth.x) ) ? (universeProperties.size.x * universeProperties.cellWidth.x - fabs(p.position.x)) : ( ( (pos.x == universeProperties.size.x - 1) && (relPos.x > halfWidth.x) ) ? p.position.x - universeProperties.size.x * universeProperties.cellWidth.x : p.position.x );

				int j = (relPos.y < -halfWidth.y) ? ( pos.y == 0 ? size.y - 1 : 0 ) : ( (relPos.y > halfWidth.y) ? ( (pos.y == universeProperties.size.y - 1) ? 0 : 2 ) : 1 );
				p.position.y = ( (pos.y == 0) && (relPos.y < -halfWidth.y) ) ? (universeProperties.size.y * universeProperties.cellWidth.y - fabs(p.position.y)) : ( ( (pos.y == universeProperties.size.y - 1) && (relPos.y > halfWidth.y) ) ? p.position.y - universeProperties.size.y * universeProperties.cellWidth.y : p.position.y );

				int k = (relPos.z < -halfWidth.z) ? ( pos.z == 0 ? size.z - 1 : 0 ) : ( (relPos.z > halfWidth.z) ? ( (pos.z == universeProperties.size.z - 1) ? 0 : 2 ) : 1 );
				p.position.z = ( (pos.z == 0) && (relPos.z < -halfWidth.z) ) ? (universeProperties.size.z * universeProperties.cellWidth.z - fabs(p.position.z)) : ( ( (pos.z == universeProperties.size.z - 1) && (relPos.z > halfWidth.z) ) ? p.position.z - universeProperties.size.z * universeProperties.cellWidth.z : p.position.z );

				// send to neighbor cell
				auto target = neighbors[{i,j,k}];
				// to place particles that exit the domain into the proper buffers
				if ( ((relPos.x < -halfWidth.x) && (pos.x == 0)) || ((relPos.y < -halfWidth.y) && (pos.y == 0)) || ((relPos.z < -halfWidth.z) && (pos.z == 0)) || ((relPos.x > halfWidth.x) && (pos.x == universeProperties.size.x - 1)) || ((relPos.y > halfWidth.y) && (pos.y == universeProperties.size.y - 1)) || ((relPos.z > halfWidth.z) && (pos.z == universeProperties.size.z - 1)) ) {
					i = ( (i == 0) || (i == size.x - 1) ) ? i : ( pos.x + (i - 1) ) * 3 + 1;
					j = ( (j == 0) || (j == size.y - 1) ) ? j : ( pos.y + (j - 1) ) * 3 + 1;
					k = ( (k == 0) || (k == size.z - 1) ) ? k : ( pos.z + (k - 1) ) * 3 + 1;
					target = &transfers[{i,j,k}];
				}
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
 	* This function verifies whether the position of the particle is within the current cell boundaries
 	*
	* @param universeProperties the properties of this universe
 	* @param cell the current cell
	* @param pos the coordinates of this cell in the grid
 	*/
	bool VerifyCorrectParticlesPositionInCell(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos) {
		int incorrectlyPlacedParticles = 0;

		for(auto& p : cell.particles) {
			// compute relative position
			Vector3<double> relPos = p.position - getCenterOfCell(pos, universeProperties);
			auto halfWidth = universeProperties.cellWidth / 2;
			if ((fabs(relPos.x) > halfWidth.x) || (fabs(relPos.y) > halfWidth.y) || (fabs(relPos.z) > halfWidth.z)) {
				++incorrectlyPlacedParticles;
			}
		}
		
		if (incorrectlyPlacedParticles) {
			std::cout << "There are " << incorrectlyPlacedParticles << " incorrectly placed particles in a cell at the position " << pos << "\n";
			return false;
		}

		return true;
	}  

	/**
	* This function imports all particles that are directed towards the specified cell
	* into the cell from the provided transfer buffers.
	*
	* @param universeProperties the properties of this universe
	* @param cell the cell to import particles into
	* @param pos the coordinates of this cell in the grid
	* @param transfers a grid of buffers to import particles from
	*/
	void importParticles(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, allscale::api::user::data::Grid<std::vector<Particle>,3>& transfers) {

		assert_true(pos.dominatedBy(universeProperties.size)) << "Position " << pos << " is outside universe of size " << universeProperties.size;

		// import particles sent to this cell
		utils::Coordinate<3> size = transfers.size();
		utils::Coordinate<3> centerIndex = pos * 3 + utils::Coordinate<3>{1,1,1};
		for(int i = 0; i<3; i++) {
			for(int j = 0; j<3; j++) {
				for(int k = 0; k<3; k++) {
					// iterates through all buffers attached to the cell at the given position
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

		// verify correct placement of the particles
		// TODO: this potential should only be used in tests due to the performance reasons
		VerifyCorrectParticlesPositionInCell(universeProperties, cell, pos);
	}

} // end namespace ipic3d
