#pragma once

#include <vector>
#include <random>

#include "allscale/api/core/io.h"
#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/operator/pfor.h"
#include "allscale/api/user/operator/ops.h"
#include "allscale/utils/static_grid.h"

#include "ipic3d/app/vector.h"
#include "ipic3d/app/particle.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/universe_properties.h"
#include "ipic3d/app/utils/points.h"
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


	Cells initCells(const Parameters& params, const InitProperties& initProperties, const UniverseProperties& properties) {

		// -- initialize the grid of cells --

		// the 3-D grid of cells
		Cells cells(properties.size);

		// -- initialize the state of each individual cell --

		// compute number of particles to be added for the uniform distribution
		unsigned particlesPerCell = initProperties.particlesPerCell[0].x * initProperties.particlesPerCell[0].y * initProperties.particlesPerCell[0].z;

		const utils::Coordinate<3> zero = 0;							// a zero constant (coordinate [0,0,0])

		// pre-compute values for computing q
		double q_factor = params.qom[0] / fabs(params.qom[0]);
		q_factor = q_factor * (properties.cellWidth.x * properties.cellWidth.y * properties.cellWidth.z) / particlesPerCell;
		
		// TODO: return this as a treeture
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {

			Cell& cell = cells[pos];

			// -- add particles --

			// add the requested number of parameters
			std::minstd_rand randGenerator((unsigned)(pos[0] * 10000 + pos[1] * 100 + pos[2]));
			const double randMax = std::minstd_rand::max();
			// TODO: can we use pfor here?
			// Maxellian random velocity and uniform spatial distribution
			// To note: this is for one spacy
			for (unsigned i = 0; i < particlesPerCell; i++) {
				Particle p;

				Vector3<double> randVals = {(double)randGenerator() / randMax, (double)randGenerator() / randMax, (double)randGenerator() / randMax};
				// initialize particle's position
				p.position = getCenterOfCell(pos, properties) + allscale::utils::elementwiseProduct(randVals, properties.cellWidth) - properties.cellWidth/2;

				// initialize particle's speed
				auto theta = 2.0 * M_PI * randVals;
				Vector3<double> prob;
				prob[0] = sqrt( -2.0 * log( 1.0 - 0.999999 * randVals[0] ) );
				prob[1] = sqrt( -2.0 * log( 1.0 - 0.999999 * randVals[1] ) );
				prob[2] = sqrt( -2.0 * log( 1.0 - 0.999999 * randVals[2] ) );

				p.velocity[0] = params.u0[0] + params.uth[0] * ( prob[0] * cos(theta[0]) );
				p.velocity[1] = params.v0[0] + params.vth[0] * ( prob[1] * sin(theta[1]) );
				p.velocity[2] = params.w0[0] + params.wth[0] * ( prob[2] * cos(theta[2]) );
				
				p.qom = params.qom[0];
				p.q = q_factor * params.rhoInit[0];  

				cell.particles.push_back(p);
			}

		});

		// return the initialized cells
		return cells;
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
		Vector3<double> Js[2][2][2] = { { { Vector3<double>(0.0) } } };

		// aggregate charge density from particles
		// TODO: use pfor here
		const auto cellOrigin = getOriginOfCell(pos, universeProperties);
		for(const auto& p : cell.particles) {
			// get the fractional distance of the particle from the cell origin
			const auto relPos = allscale::utils::elementwiseDivision((p.position - cellOrigin), (universeProperties.cellWidth));

			// computation of J also includes weights from the particles as for E
			for(int i = 0; i < 2; ++i) {
				for(int j = 0; j < 2; ++j) {
					for(int k = 0; k < 2; ++k) {
						//utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});
						//auto cornerPos = getOriginOfCell(cur, universeProperties); 
						//auto fac = (i == 0 ? (1 - cornerPos.x) : cornerPos.x) * (j == 0 ? (1 - cornerPos.y) : cornerPos.y) * (k == 0 ? (1 - cornerPos.z) : cornerPos.z);
				    	auto fac = (i == 0 ? (1 - relPos.x) : relPos.x) * (j == 0 ? (1 - relPos.y) : relPos.y) * (k == 0 ? (1 - relPos.z) : relPos.z);
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

		const auto cellOrigin = getOriginOfCell(pos, universeProperties);

        double vol = universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;

		// update particles
		allscale::api::user::pfor(cell.particles, [&](Particle& p){
			// Docu: https://www.particleincell.com/2011/vxb-rotation/
			// Code: https://www.particleincell.com/wp-content/uploads/2011/07/ParticleIntegrator.java

			// get the fractional distance of the particle from the cell origin
			const auto relPos = allscale::utils::elementwiseDivision((p.position - cellOrigin), (universeProperties.cellWidth));

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
		allscale::utils::StaticGrid<std::vector<Particle>*,3,3,3> neighbors;
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
			auto halfWidth = universeProperties.cellWidth / 2.0;
			if((fabs(relPos.x) > halfWidth.x) || (fabs(relPos.y) > halfWidth.y) || (fabs(relPos.z) > halfWidth.z)) {
				std::cout << "New particle\n";
				std::cout << pos << '\n';
				std::cout << p.position << '\n';
				std::cout << relPos << '\n';

				// compute corresponding neighbor cell
				// cover the inner cells as well as the boundary cells on positions 0 and N-1
				auto computeCell = [&](const int i) {
					return (relPos[i] < -halfWidth[i]) ? (pos[i] == 0 ? size[i] - 1 : 0) : ((relPos[i] > halfWidth[i]) ? ((pos[i] == universeProperties.size[i] - 1) ? 0 : 2) : 1);
				};

				// adjust particle's position in case it exits the domain
				auto adjustPosition = [&](const int i) {
					return p.position[i] = ((pos[i] == 0) && (relPos[i] < -halfWidth[i])) ? (universeProperties.origin[i] + universeProperties.size[i] * universeProperties.cellWidth[i] - fabs(p.position[i])) : (((pos[i] == universeProperties.size[i] - 1) && (relPos[i] > halfWidth[i])) ? p.position[i] - (universeProperties.origin[i] + universeProperties.size[i] * universeProperties.cellWidth[i]) : p.position[i]);
				};

				int i = (int)computeCell(0);
				int j = (int)computeCell(1);
				int k = (int)computeCell(2);
				std::cout << i << '\t';
				std::cout << j << '\t';
				std::cout << k << '\n';

				p.position[0] = adjustPosition(0);
				p.position[1] = adjustPosition(1);
				p.position[2] = adjustPosition(2);
//				std::cout << p.position << '\n';
//				std::cout << "Position\n";

				// send to neighbor cell
				auto target = neighbors[{i,j,k}];
				// to place particles that exit the domain into the proper buffers
				if ( ((relPos.x < -halfWidth.x) && (pos.x == 0)) || ((relPos.y < -halfWidth.y) && (pos.y == 0)) || ((relPos.z < -halfWidth.z) && (pos.z == 0)) || ((relPos.x > halfWidth.x) && (pos.x == universeProperties.size.x - 1)) || ((relPos.y > halfWidth.y) && (pos.y == universeProperties.size.y - 1)) || ((relPos.z > halfWidth.z) && (pos.z == universeProperties.size.z - 1)) ) {
					i = ( (i == 0) || (i == (int)size.x - 1) ) ? i : ( (int)pos.x + (i - 1) ) * 3 + 1;
					j = ( (j == 0) || (j == (int)size.y - 1) ) ? j : ( (int)pos.y + (j - 1) ) * 3 + 1;
					k = ( (k == 0) || (k == (int)size.z - 1) ) ? k : ( (int)pos.z + (k - 1) ) * 3 + 1;
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
			auto halfWidth = universeProperties.cellWidth / 2.0;
			if ((fabs(relPos.x) > halfWidth.x) || (fabs(relPos.y) > halfWidth.y) || (fabs(relPos.z) > halfWidth.z)) {
				++incorrectlyPlacedParticles;
				std::cout << "Error \n";
				std::cout << pos << '\n';
				std::cout << getOriginOfCell(pos, universeProperties) << '\n';
				std::cout << p.position << '\n';
                exit(13);	
            }
		}
		
		if (incorrectlyPlacedParticles) {
			std::cout << "There are " << incorrectlyPlacedParticles << " incorrectly placed particles in a cell at the position " << pos << "\n";
			incorrectlyPlacedParticles = 0;
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

	/** 
 	 * This function computes particles total kinetic energy
 	 */
	double getParticlesKineticEnergy(Cell& cell) {
		// TODO: use more convenient reduction operators once they are available in the API
		auto map = [](const Particle& p, double& res) {
			res += 0.5 * (p.q / p.qom) * allscale::utils::sumOfSquares(p.velocity);
		};

		auto reduce = [&](const double& a, const double& b) { return a + b; };
		auto init = []() { return 0.0; };

		return allscale::api::user::preduce(cell.particles, map, reduce, init);
	} 

	/** 
 	 * This function computes particles total momentum
 	 */
	double getParticlesMomentum(Cell& cell) {
		// TODO: use more convenient reduction operators once they are available in the API
		auto map = [](const Particle& p, double& res) {
			res += (p.q / p.qom) * sqrt(allscale::utils::sumOfSquares(p.velocity)); 
		};

		auto reduce = [&](const double& a, const double& b) { return a + b; };
		auto init = []() { return 0.0; };

		return allscale::api::user::preduce(cell.particles, map, reduce, init);
	}

	/**
	* This function outputs the number of particles per cell
	*/
	template<typename StreamObject>
	void outputNumberOfParticlesPerCell(const Cells& cells, StreamObject& streamObject) {
		// TODO: implement output facilities for large problems
		assert_le(cells.size(), (coordinate_type{ 32,32,32 })) << "Unable to dump data for such a large cell grid at this time";

		// output dimensions
		streamObject << cells.size() << "\n";

		// output particles per cell
		allscale::api::user::pfor(cells.size(), [&](const auto& index) {
			streamObject.atomic([&](auto& out) { 
				out << index.x << "," << index.y << "," << index.z << ":";
				out << cells[index].particles.size() << "\n"; });
		});

		streamObject << "\n";
	}

	/**
	* This function outputs the number of particles per cell using AllScale IO
	*/
	void outputNumberOfParticlesPerCell(const Cells& cells, std::string& filename) {
		auto& manager = allscale::api::core::FileIOManager::getInstance();
		auto text = manager.createEntry(filename);
		auto out = manager.openOutputStream(text);
		outputNumberOfParticlesPerCell(cells, out);
		manager.close(out);
	}

} // end namespace ipic3d
