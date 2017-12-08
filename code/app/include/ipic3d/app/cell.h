#pragma once

#include <vector>
#include <random>

#include "allscale/api/core/io.h"
#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/user/algorithm/preduce.h"
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
		unsigned totalParticlesPerCell = initProperties.particlesPerCell[0].x * initProperties.particlesPerCell[0].y * initProperties.particlesPerCell[0].z;

		const utils::Coordinate<3> zero = 0;							// a zero constant (coordinate [0,0,0])

		// pre-compute values for computing q
		double q_factor = params.qom[0] / fabs(params.qom[0]);
		q_factor = q_factor * (properties.cellWidth.x * properties.cellWidth.y * properties.cellWidth.z) / totalParticlesPerCell;
		const double fourPI = 16.0 * atan(1.0);
		q_factor = q_factor * params.rhoInit[0] / fourPI;

		auto particlesPerCell = initProperties.particlesPerCell[0];
	
		// TODO: return this as a treeture
		allscale::api::user::algorithm::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {

			Cell& cell = cells[pos];
			auto cellOrigin = getOriginOfCell(pos, properties);

			// add the requested number of parameters
			std::minstd_rand randGenerator((unsigned)(pos[0] * 10000 + pos[1] * 100 + pos[2]));
			const double randMax = std::minstd_rand::max();

			// -- add particles --
			// TODO: we plan to can use bags to store particle which would allow us to parallelize this for loop
			// Maxellian random velocity and uniform spatial distribution
			// To note: this is for one specie
			for (unsigned i = 0; i < particlesPerCell.x; i++) {
				for (unsigned j = 0; j < particlesPerCell.y; j++) {
					for (unsigned k = 0; k < particlesPerCell.z; k++) {
						Particle p;

						// initialize particle's position
						p.position.x = (i + 0.5) * (properties.cellWidth.x / particlesPerCell.x) + cellOrigin.x; 
						p.position.y = (j + 0.5) * (properties.cellWidth.y / particlesPerCell.y) + cellOrigin.y; 
						p.position.z = (k + 0.5) * (properties.cellWidth.z / particlesPerCell.z) + cellOrigin.z; 

						// initialize particle's velocity
						double prob0, prob1;
						double theta0, theta1;

						double harvest = (double)randGenerator() / randMax;
						prob0 = sqrt( -2.0 * log( 1.0 - 0.999999 * harvest ) );
						harvest = (double)randGenerator() / randMax;
						theta0 = 2.0 * M_PI * harvest;

						harvest = (double)randGenerator() / randMax;
						prob1 = sqrt( -2.0 * log( 1.0 - 0.999999 * harvest ) );
						harvest = (double)randGenerator() / randMax;
						theta1 = 2.0 * M_PI * harvest;

						p.velocity.x = params.u0[0] + params.uth[0] * ( prob0 * cos(theta0) );
						p.velocity.y = params.v0[0] + params.vth[0] * ( prob0 * sin(theta0) );
						p.velocity.z = params.w0[0] + params.wth[0] * ( prob1 * cos(theta1) );
						
						p.qom = params.qom[0];
						p.q = q_factor;  

						cell.particles.push_back(p);
					}
				}
			}
	
			// print particles position and velocity
			if (0) {
				for(const auto& p : cell.particles) {
					std::cout << p.position.x << " ";
					std::cout << p.position.y << " ";
					std::cout << p.position.z << " ";
					std::cout << p.velocity.x << " ";
					std::cout << p.velocity.y << " ";
					std::cout << p.velocity.z << "\n";
				}
			
			}

		});

		// return the initialized cells
		return cells;
	}

	/**
	* This function projects the effect of the particles contained in the specified cell
	* to the given density contributions defined on the grid
	*
	* @param universeProperties the properties of this universe
	* @param cell the cell whose particles are considered in the density contributions computation
	* @param pos the coordinates of this cell in the grid
	* @param contributions the density contributions output
	*/
	void projectToDensityField(const UniverseProperties& universeProperties, const Cells& cells, const utils::Coordinate<3>& pos, CurrentDensity& density) {

		auto Js = Vector3<double>(0.0);
		auto size = universeProperties.size;

		// aggregate charge density from particles
		for(int i=-1; i<1; i++) {
			for(int j=-1; j<1; j++) {
				for(int k=-1; k<1; k++) {
					utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});

					// adjust particle's position in case it exits the domain
					auto adjustPosition = [&](const int i) {
						return cur[i] = ( (cur[i] < 0) ? size[i] - 1 : ((cur[i] >= size[i]) ? 0 : cur[i]) );
					};
				
					cur[0] = adjustPosition(0);
					cur[1] = adjustPosition(1);
					cur[2] = adjustPosition(2);

					const auto cellOrigin = getOriginOfCell(cur, universeProperties);
					for(const auto& p : cells[cur].particles) {
						// get the fractional distance of the particle from the cell origin
						const auto relPos = allscale::utils::elementwiseDivision((p.position - cellOrigin), (universeProperties.cellWidth));

						// computation of J also includes weights from the particles as for E
						// despite the fact that we are working right now with multiple cells, so the position of J would be different
						// 	the formula still works well as it captures position of J in each of those cells. 
				    	auto fac = (i == 0 ? (1 - relPos.x) : relPos.x) * (j == 0 ? (1 - relPos.y) : relPos.y) * (k == 0 ? (1 - relPos.z) : relPos.z);
						Js += p.q * p.velocity * fac;
					}
				}
			}
		}

		double vol = universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;
		density[pos].J = (Js / vol) / 8.0;
	}

	/**
	* This function computes a trilinear interpolation for a given target position within a rectangular box spanned by 8 corner points
	*   This is interpolation of fields to particles
	* Math: http://paulbourke.net/miscellaneous/interpolation/
	* TODO: Move this to some math utilities header?
			etd::cout << v_plus << '\n';
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
	void moveParticles(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, const Field& field) {

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
		allscale::api::user::algorithm::pfor(cell.particles, [&](Particle& p){
			// Docu: https://www.particleincell.com/2011/vxb-rotation/
			// Code: https://www.particleincell.com/wp-content/uploads/2011/07/ParticleIntegrator.java

			// get the fractional distance of the particle from the cell origin
			const auto relPos = allscale::utils::elementwiseDivision((p.position - cellOrigin), (universeProperties.cellWidth));

			// interpolate
			auto E = trilinearInterpolationF2P(Es, relPos, vol);
			auto B = trilinearInterpolationF2P(Bs, relPos, vol);

			// update velocity
			p.updateVelocity(E, B, universeProperties.dt);

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
					// wraps indices both for overflow and underflow (-1 is mapped to size - 1)
					// implementing periodic boundary conditions
					auto wrap = [](const auto value, const auto delta, const auto max) {
						if (delta >= 0)
							return (value + delta) % max;
						else
							return ((value + delta) + max * (-delta)) % max;
					};

					const auto offset = utils::Coordinate<3>{ i - 1,j - 1,k - 1 } *2;
					auto targetPos = centerIndex + offset;
					targetPos[0] = wrap(centerIndex[0], offset[0], size[0]);
					targetPos[1] = wrap(centerIndex[1], offset[1], size[1]);
					targetPos[2] = wrap(centerIndex[2], offset[2], size[2]);

					neighbors[{i,j,k}] = &transfers[targetPos];
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

				// adjust particle's position in case it exits the domain
				auto adjustPosition = [&](const int i) {
					return p.position[i] = ((pos[i] == 0) && (relPos[i] < -halfWidth[i])) ? (universeProperties.size[i] * universeProperties.cellWidth[i] + p.position[i]) : (((pos[i] == universeProperties.size[i] - 1) && (relPos[i] > halfWidth[i])) ? p.position[i] - universeProperties.size[i] * universeProperties.cellWidth[i] : p.position[i]);
				};

				int i = (relPos.x < -halfWidth.x) ? 0 : ((relPos.x > halfWidth.x) ? 2 : 1);
				int j = (relPos.y < -halfWidth.y) ? 0 : ((relPos.y > halfWidth.y) ? 2 : 1);
				int k = (relPos.z < -halfWidth.z) ? 0 : ((relPos.z > halfWidth.z) ? 2 : 1);

				p.position[0] = adjustPosition(0);
				p.position[1] = adjustPosition(1);
				p.position[2] = adjustPosition(2);

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
 	* This function verifies whether the position of the particle is within the current cell boundaries
 	*
	* @param universeProperties the properties of this universe
 	* @param cell the current cell
	* @param pos the coordinates of this cell in the grid
 	*/
	bool verifyCorrectParticlesPositionInCell(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos) {
		int incorrectlyPlacedParticles = 0;

		for(auto& p : cell.particles) {
			// compute relative position
			Vector3<double> relPos = p.position - getCenterOfCell(pos, universeProperties);
			auto halfWidth = universeProperties.cellWidth / 2.0;
			if ((fabs(relPos.x) > halfWidth.x) || (fabs(relPos.y) > halfWidth.y) || (fabs(relPos.z) > halfWidth.z)) {
				++incorrectlyPlacedParticles;
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
		assert_true(verifyCorrectParticlesPositionInCell(universeProperties, cell, pos));
	}

	/** 
 	 * This function computes particles total kinetic energy in a cell
 	 */
	double getParticlesKineticEnergy(const Cell& cell) {
		auto map = [](const Particle& p, double& res) {
			res += 0.5 * (p.q / p.qom) * allscale::utils::sumOfSquares(p.velocity);
		};

		auto reduce = [&](const double& a, const double& b) { return a + b; };
		auto init = []() { return 0.0; };

		return allscale::api::user::algorithm::preduce(cell.particles, map, reduce, init).get();
	} 

	/** 
 	 * This function computes particles total momentum in a cell
 	 */
	double getParticlesMomentum(const Cell& cell) {
		auto map = [](const Particle& p, double& res) {
			res += (p.q / p.qom) * sqrt(allscale::utils::sumOfSquares(p.velocity)); 
		};

		auto reduce = [&](const double& a, const double& b) { return a + b; };
		auto init = []() { return 0.0; };

		return allscale::api::user::algorithm::preduce(cell.particles, map, reduce, init).get();
	}

	/**
 	 * This functions computes total particles energy.
 	 * 	This is the second reduction among values computed on the cells level
 	 * 	The first reduction is the computation of particles energies on the cells level
 	 */
	template<class T>
	double getTotalParticlesEnergy(const Cells& cells, T func){
		auto map = [&](const coordinate_type& index, double& res) {
			res += func(cells[index]);
		};

		auto reduce = [&](const double& a, const double& b) { return a + b; };
		auto init = []() { return 0.0; };

		return allscale::api::user::algorithm::preduce(coordinate_type(0), cells.size(), map, reduce, init).get();
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
		allscale::api::user::algorithm::pfor(cells.size(), [&](const auto& index) {
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
