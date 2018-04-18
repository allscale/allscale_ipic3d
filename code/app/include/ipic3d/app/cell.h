#pragma once

#include <vector>
#include <random>

#include "allscale/api/core/io.h"
#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/user/algorithm/async.h"
#include "allscale/api/user/algorithm/preduce.h"
#include "allscale/utils/static_grid.h"

#include "ipic3d/app/field.h"
#include "ipic3d/app/parameters.h"
#include "ipic3d/app/particle.h"
#include "ipic3d/app/transfer_buffer.h"
#include "ipic3d/app/universe_properties.h"
#include "ipic3d/app/utils/points.h"
#include "ipic3d/app/vector.h"
#include "ipic3d/app/ziggurat_normal_distribution.h"

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


	/**
	 * Tests whether a given particle is to be maintained by a cell of the given position.
	 * @param universeProperties the properties of this universe
	 * @param pos the coordinates of this cell in the grid
	 * @param p the particle to be tested
	 * @return true if the particle should be stored in the designated cell, false otherwise
	 */
	bool isInside(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, const Particle& p) {
		// compute relative position
		Vector3<double> relPos = p.position - getCenterOfCell(pos, universeProperties);

		// get cell width
		auto halfWidth = universeProperties.cellWidth / 2.0;

		// TODO: sort out which face belongs to which cell

		// check relative position
		return ((fabs(relPos.x) <= halfWidth.x) && (fabs(relPos.y) <= halfWidth.y) && (fabs(relPos.z) <= halfWidth.z));
	}

	/**
	 * Tests whether a given particle is within the particle universe.
	 * @param universeProperties the properties of this universe
	 * @param p the particle to be tested
	 * @return true if the particle is inside, false otherwise
	 */
	bool isInsideUniverse(const UniverseProperties& universeProperties, const Particle& p) {
		Vector3<double> zero = 0;
		Vector3<double> size {
			universeProperties.size.x * universeProperties.cellWidth.x,
			universeProperties.size.y * universeProperties.cellWidth.y,
			universeProperties.size.z * universeProperties.cellWidth.z
		};
		return zero.dominatedBy(p.position) && p.position.strictlyDominatedBy(size);
	}


	// count the number of particles in all cells
	std::uint64_t countParticlesInDomain(const Cells& cells) {

		auto fold = [&](const coordinate_type& index, std::uint64_t& res) {
			res += (std::uint64_t)cells[index].particles.size();
		};

		auto reduce = [&](const std::uint64_t& a, const std::uint64_t& b) { return a + b; };
		auto init = []()->std::uint64_t { return 0; };

		coordinate_type zero(0);
		coordinate_type full(cells.size());

		return allscale::api::user::algorithm::preduce(zero, full, fold, reduce, init).get();
	}

	namespace distribution {

		namespace species {

			// a generator for electrons
			struct electron {

				Particle operator()() const {
					Particle p;
					p.q = -1.0;
					p.qom = -25.0;
					return p;
				}

			};

			// a generator for protons
			struct proton {

				Particle operator()() const {
					Particle p;
					p.q = 1.0;
					p.qom = 1.0;
					return p;
				}

			};

		}

		namespace vector {

			// a generator for uniformly distributed vector3 instances
			class uniform {

				std::uniform_real_distribution<> x;
				std::uniform_real_distribution<> y;
				std::uniform_real_distribution<> z;

				std::minstd_rand randGen;

			public:

				uniform(const Vector3<double>& min, const Vector3<double>& max, std::uint32_t seed)
					: x(min.x,max.x), y(min.y,max.y), z(min.z,max.z), randGen(seed) {}

				Vector3<double> operator()() {
					return { x(randGen), y(randGen), z(randGen) };
				}

			};

			// a generator for normal distributed vector3 instances
			class normal {

				Vector3<double> mean;
				Vector3<double> stddev;

				ziggurat_normal_distribution rand;

			public:

				normal(const Vector3<double>& mean, const Vector3<double>& stddev, std::uint32_t seed)
					: mean(mean), stddev(stddev), rand(seed) {}

				Vector3<double> operator()() {
					return {
						mean.x + stddev.x * rand(),
						mean.y + stddev.y * rand(),
						mean.z + stddev.z * rand()
					};
				}

			};

		}


		template<typename PositionDist, typename VelocityDist, typename SpeciesDist>
		class generic_particle_generator {

			SpeciesDist speciesGen;
			PositionDist posGen;
			VelocityDist velGen;

		public:

			generic_particle_generator(const PositionDist& posGen, const VelocityDist& velGen, const SpeciesDist& speciesGen)
				: speciesGen(speciesGen), posGen(posGen), velGen(velGen) {}

			Particle operator()() {
				Particle p = speciesGen();
				p.position = posGen();
				p.velocity = velGen();
				return p;
			}

		};


		template<typename SpeciesGen = species::electron>
		struct uniform : public generic_particle_generator<vector::uniform,vector::uniform,SpeciesGen> {

			using super = generic_particle_generator<vector::uniform,vector::uniform,SpeciesGen>;

			uniform(
					const Vector3<double>& minPos,
					const Vector3<double>& maxPos,
					const Vector3<double>& minVel,
					const Vector3<double>& maxVel,
					std::uint32_t seed = 0
			) : super({minPos,maxPos,seed+1},{minVel,maxVel,seed+2},SpeciesGen()) {}

			uniform(
					const SpeciesGen& speciesGen,
					const Vector3<double>& minPos,
					const Vector3<double>& maxPos,
					const Vector3<double>& minVel,
					const Vector3<double>& maxVel,
					std::uint32_t seed = 0
			) : super({minPos,maxPos,seed+1},{minVel,maxVel,seed+2},speciesGen) {}

		};

		template<typename SpeciesGen = species::electron>
		struct normal : public generic_particle_generator<vector::normal,vector::uniform,SpeciesGen> {

			using super = generic_particle_generator<vector::normal,vector::uniform,SpeciesGen>;

			normal(
					const Vector3<double>& center,
					const Vector3<double>& stddev,
					const Vector3<double>& minVel,
					const Vector3<double>& maxVel,
					std::uint32_t seed = 0
			) : super({center,stddev,seed+1},{minVel,maxVel,seed+2},SpeciesGen()) {}

			normal(
					const SpeciesGen& speciesGen,
					const Vector3<double>& center,
					const Vector3<double>& stddev,
					const Vector3<double>& minVel,
					const Vector3<double>& maxVel,
					std::uint32_t seed = 0
			) : super({center,stddev,seed+1},{minVel,maxVel,seed+2},speciesGen) {}

		};

		template<typename Distribution>
		class spherical {

			Distribution dist;

			Vector3<double> center;
			double radius;

		public:

			spherical(const Distribution& dist, const Vector3<double>& center, double radius)
				: dist(dist), center(center), radius(radius) {
				assert_gt(radius,0);
			}

			Particle operator()() {
				for(int i=0; i<100000; i++) {
					Particle p = dist();
					if (norm(p.position - center) <= radius) {
						return p;
					}
				}
				assert_fail() << "Probably empty distribution!";
				return {};
			}

		};

		template<typename Distribution>
		spherical<Distribution> make_spherical(const Distribution& dist, const Vector3<double>& center = 0, double radius = 1) {
			return { dist, center, radius };
		}

	}


	template<typename Distribution>
	Cells initCells(const UniverseProperties& properties, std::uint64_t numParticles, const Distribution& dist) {
		using allscale::api::user::algorithm::pfor;

		// the 3-D grid of cells
		Cells cells(properties.size);

		// get a private copy of the distribution generator
		auto next = dist;

		// just some info about the progress
		std::cout << "Sorting in particles ...\n";

		// insert particles in patches (of 4MB each)
		const std::uint64_t patch_size = (1<<22)/sizeof(Particle);
		for(std::uint64_t p=0; p<numParticles; p+=patch_size) {

			// generate list of particles
			auto current_patch = std::min(patch_size,numParticles-p);
			std::vector<Particle> particles(current_patch);
			for(std::uint64_t i=0; i<current_patch; ++i) {
				auto cur = next();
				while(!isInsideUniverse(properties,cur)) cur = next();
				particles[i] = cur;
			}

			std::cout << "Submitting particles " << p << " - " << (p+current_patch) << " ... \n";

			// distribute particles randomly
			pfor(properties.size, [=,&cells](const auto& pos) {

				// get targeted cell
				auto& cell = cells[pos];

				// get cell corners
				auto& width = properties.cellWidth;
				Vector3<double> low { width.x * pos.x, width.y * pos.y, width.z * pos.z };
				Vector3<double> hig = low + width;

				// filter out local particles
				for(const auto& cur : particles) {
					if (low.dominatedBy(cur.position) && cur.position.strictlyDominatedBy(hig)) {
						cell.particles.push_back(cur);
					}
				}

			});

		}

		// done
		return cells;
	}

	Cells initCells(const UniverseProperties& properties, std::uint64_t numParticles, const distribution::uniform<>&) {
		using allscale::api::user::algorithm::pfor;

		// the 3-D grid of cells
		auto gridSize = properties.size;
		Cells cells(gridSize);

		// just some info about the progress
		std::cout << "Sorting in uniformly distributed particles ...\n";

		// compute particles per cell
		auto numCells = properties.size.x * properties.size.y * properties.size.z;
		std::uint64_t particlesPerCell = numParticles / numCells;
		std::uint64_t remaining = numParticles % numCells;

		std::cout << "  particles / cell: " << particlesPerCell << " (+1)\n";

		// initialize each cell in parallel
		pfor(properties.size, [=,&cells](const auto& pos) {

			// get targeted cell
			auto& cell = cells[pos];

			// get cell corners
			auto& width = properties.cellWidth;
			Vector3<double> low { width.x * pos.x, width.y * pos.y, width.z * pos.z };
			Vector3<double> hig = low + width;

			// create a uniform distribution
			auto seed = ((pos.x * 1023) + pos.y) * 1023 + pos.z;

			// TODO: speeds are hard-coded, actual passed distribution is ignored
			distribution::uniform<> next(
					low,hig, // within this box
					// speeds are constant
					Vector3<double> { -0.2, -0.2, -0.2},
					Vector3<double> { +0.2, +0.2, +0.2},
					seed
			);

			// get number of particles to be generated in this cell
			auto numParticles = particlesPerCell;

			// correct for remaining particles (to be evenly balance)
			std::uint64_t linPos = (pos.x * gridSize.y + pos.y) * gridSize.z + pos.z;
			if (linPos < remaining) {
				numParticles += 1;
			}

			// generate particles
			for(std::uint64_t i=0; i<numParticles; i++) {
				cell.particles.push_back(next());
			}

		});

		// done
		return cells;
	}


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
		allscale::api::user::algorithm::pfor(zero, properties.size, [=,&cells](const utils::Coordinate<3>& pos) {

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
	
//			// print particles position and velocity
//			if (0) {
//				for(const auto& p : cell.particles) {
//					std::cout << p.position.x << " ";
//					std::cout << p.position.y << " ";
//					std::cout << p.position.z << " ";
//					std::cout << p.velocity.x << " ";
//					std::cout << p.velocity.y << " ";
//					std::cout << p.velocity.z << "\n";
//				}
//
//			}

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
					utils::Coordinate<3> cur({pos[0]+i+1,pos[1]+j+1,pos[2]+k+1});
					Es[i][j][k] = field[cur].E;
					Bs[i][j][k] = field[cur].B;
				}
			}
		}

		const auto cellOrigin = getOriginOfCell(pos, universeProperties);

		double vol = universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;

		// update particles
//		allscale::api::user::algorithm::pfor(cell.particles, [&](Particle& p){
		for(std::size_t i=0; i<cell.particles.size(); ++i) {
			Particle& p = cell.particles[i];
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
		}
//		});

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
	void exportParticles(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, TransferBuffers& transfers) {

		assert_true(pos.dominatedBy(universeProperties.size)) << "Position " << pos << " is outside universe of size " << universeProperties.size;

		// -- migrate particles to other cells if boundaries are crossed --

		// create buffer of remaining particles
		std::vector<Particle> remaining;
		remaining.reserve(cell.particles.size());

		{

			auto size = universeProperties.size;

			std::vector<Particle>* neighbors[3][3][3];

			// NOTE: due to an unimplemented feature in the analysis, this loop needs to be unrolled (work in progress)

//			for(int i = 0; i<3; i++) {
//				for(int j = 0; j<3; j++) {
//					for(int k = 0; k<3; k++) {
//
//						// get neighbor cell in this direction (including wrap-around)
//						auto neighbor = (pos + utils::Coordinate<3>{ i-1, j-1, k-1 } + size) % size;
//
//						// index this buffer to the neighboring cell
//						auto& buffer = transfers.getBuffer(neighbor,TransferDirection{ 2-i, 2-j, 2-k });
//						neighbors[i][j][k] = &buffer;
//
//						// clear buffer from particles moved in the last iteration
//						buffer.clear();
//					}
//				}
//			}

			// -- unroll begin --

			neighbors[0][0][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, -1, -1 } + size) % size,TransferDirection{ 2, 2, 2 });
			neighbors[0][0][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, -1,  0 } + size) % size,TransferDirection{ 2, 2, 1 });
			neighbors[0][0][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, -1,  1 } + size) % size,TransferDirection{ 2, 2, 0 });

			neighbors[0][1][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1,  0, -1 } + size) % size,TransferDirection{ 2, 1, 2 });
			neighbors[0][1][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1,  0,  0 } + size) % size,TransferDirection{ 2, 1, 1 });
			neighbors[0][1][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1,  0,  1 } + size) % size,TransferDirection{ 2, 1, 0 });

			neighbors[0][2][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1,  1, -1 } + size) % size,TransferDirection{ 2, 0, 2 });
			neighbors[0][2][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1,  1,  0 } + size) % size,TransferDirection{ 2, 0, 1 });
			neighbors[0][2][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1,  1,  1 } + size) % size,TransferDirection{ 2, 0, 0 });


			neighbors[1][0][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, -1, -1 } + size) % size,TransferDirection{ 1, 2, 2 });
			neighbors[1][0][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, -1,  0 } + size) % size,TransferDirection{ 1, 2, 1 });
			neighbors[1][0][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, -1,  1 } + size) % size,TransferDirection{ 1, 2, 0 });

			neighbors[1][1][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0,  0, -1 } + size) % size,TransferDirection{ 1, 1, 2 });
			neighbors[1][1][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0,  0,  0 } + size) % size,TransferDirection{ 1, 1, 1 });
			neighbors[1][1][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0,  0,  1 } + size) % size,TransferDirection{ 1, 1, 0 });

			neighbors[1][2][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0,  1, -1 } + size) % size,TransferDirection{ 1, 0, 2 });
			neighbors[1][2][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0,  1,  0 } + size) % size,TransferDirection{ 1, 0, 1 });
			neighbors[1][2][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0,  1,  1 } + size) % size,TransferDirection{ 1, 0, 0 });


			neighbors[2][0][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, -1, -1 } + size) % size,TransferDirection{ 0, 2, 2 });
			neighbors[2][0][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, -1,  0 } + size) % size,TransferDirection{ 0, 2, 1 });
			neighbors[2][0][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, -1,  1 } + size) % size,TransferDirection{ 0, 2, 0 });

			neighbors[2][1][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1,  0, -1 } + size) % size,TransferDirection{ 0, 1, 2 });
			neighbors[2][1][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1,  0,  0 } + size) % size,TransferDirection{ 0, 1, 1 });
			neighbors[2][1][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1,  0,  1 } + size) % size,TransferDirection{ 0, 1, 0 });

			neighbors[2][2][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1,  1, -1 } + size) % size,TransferDirection{ 0, 0, 2 });
			neighbors[2][2][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1,  1,  0 } + size) % size,TransferDirection{ 0, 0, 1 });
			neighbors[2][2][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1,  1,  1 } + size) % size,TransferDirection{ 0, 0, 0 });

			neighbors[0][0][0]->clear();
			neighbors[0][0][1]->clear();
			neighbors[0][0][2]->clear();

			neighbors[0][1][0]->clear();
			neighbors[0][1][1]->clear();
			neighbors[0][1][2]->clear();

			neighbors[0][2][0]->clear();
			neighbors[0][2][1]->clear();
			neighbors[0][2][2]->clear();


			neighbors[1][0][0]->clear();
			neighbors[1][0][1]->clear();
			neighbors[1][0][2]->clear();

			neighbors[1][1][0]->clear();
			neighbors[1][1][1]->clear();
			neighbors[1][1][2]->clear();

			neighbors[1][2][0]->clear();
			neighbors[1][2][1]->clear();
			neighbors[1][2][2]->clear();


			neighbors[2][0][0]->clear();
			neighbors[2][0][1]->clear();
			neighbors[2][0][2]->clear();

			neighbors[2][1][0]->clear();
			neighbors[2][1][1]->clear();
			neighbors[2][1][2]->clear();

			neighbors[2][2][0]->clear();
			neighbors[2][2][1]->clear();
			neighbors[2][2][2]->clear();

			// -- unroll end --

			// sort out particles
			std::vector<std::vector<Particle>*> targets(cell.particles.size());
//			allscale::api::user::algorithm::pfor(0ul,cell.particles.size(),[&](std::size_t index){
			for(std::size_t index=0; index<cell.particles.size(); ++index) {

				// get the current particle
				auto& p = cell.particles[index];

				// compute relative position
				Vector3<double> relPos = p.position - getCenterOfCell(pos, universeProperties);
				auto halfWidth = universeProperties.cellWidth / 2.0;

				// if required, "reflect" particle's position and mark that velocity vector should be inverted
				bool invertVelocity = false;
				auto adjustPosition = [&](const int i) {
					if((pos[i] == 0) && (relPos[i] < -halfWidth[i])) {
						invertVelocity = true;
						return p.position[i] + halfWidth[i];
					} else if((pos[i] == universeProperties.size[i] - 1) && (relPos[i] > halfWidth[i])) {
						invertVelocity = true;
						return p.position[i] - halfWidth[i];
					}
					return p.position[i];
				};

				p.position[0] = adjustPosition(0);
				p.position[1] = adjustPosition(1);
				p.position[2] = adjustPosition(2);

				if(invertVelocity) {
					p.velocity *= (-1);
				}

				// recompute potentially new relative position
				relPos = p.position - getCenterOfCell(pos, universeProperties);

				// send particle to neighboring cell if required
				if((fabs(relPos.x) > halfWidth.x) || (fabs(relPos.y) > halfWidth.y) || (fabs(relPos.z) > halfWidth.z)) {

					int i = (relPos.x < -halfWidth.x) ? 0 : ((relPos.x > halfWidth.x) ? 2 : 1);
					int j = (relPos.y < -halfWidth.y) ? 0 : ((relPos.y > halfWidth.y) ? 2 : 1);
					int k = (relPos.z < -halfWidth.z) ? 0 : ((relPos.z > halfWidth.z) ? 2 : 1);

					targets[index] = neighbors[i][j][k];
				} else {
					// keep particle
					targets[index] = &remaining;
				}
//			});
			}

			// actually transfer particles
			for(std::size_t i = 0; i<cell.particles.size(); ++i) {
				targets[i]->push_back(cell.particles[i]);
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
			if (!isInside(universeProperties,pos,p)) {
				++incorrectlyPlacedParticles;
			}
		}
		
		if (incorrectlyPlacedParticles) {
			std::cerr << "There are " << incorrectlyPlacedParticles << " incorrectly placed particles in a cell at the position " << pos << "\n";
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
	void importParticles(const UniverseProperties& universeProperties, Cell& cell, const utils::Coordinate<3>& pos, TransferBuffers& transfers) {

		assert_true(pos.dominatedBy(universeProperties.size)) << "Position " << pos << " is outside universe of size " << universeProperties.size;

		// import particles sent to this cell

//		for(int i = 0; i<3; i++) {
//			for(int j = 0; j<3; j++) {
//				for(int k = 0; k<3; k++) {
//
//					// skip the center (not relevant)
////					if (i == 1 && j == 1 && k == 1) continue;
//
//					// obtain transfer buffer
//					auto& in = transfers.getBuffer(pos,TransferDirection(i,j,k));
//
//					// import particles
//					cell.particles.insert(cell.particles.end(), in.begin(), in.end());
//				}
//			}
//		}


		// NOTE: due to an unimplemented feature in the analysis, this loop needs to be unrolled (work in progress)

		auto import = [&](const auto& in) {
			if (in.empty()) return;
			auto& cur = cell.particles;
			auto oldSize = cur.size();
			cur.resize(oldSize + in.size());
			std::memcpy(&cur[oldSize],&in[0],sizeof(Particle) * in.size());
		};

		// along all 26 directions (center is not relevant)
		import(transfers.getBuffer(pos,TransferDirection(0,0,0)));
		import(transfers.getBuffer(pos,TransferDirection(0,0,1)));
		import(transfers.getBuffer(pos,TransferDirection(0,0,2)));

		import(transfers.getBuffer(pos,TransferDirection(0,1,0)));
		import(transfers.getBuffer(pos,TransferDirection(0,1,1)));
		import(transfers.getBuffer(pos,TransferDirection(0,1,2)));

		import(transfers.getBuffer(pos,TransferDirection(0,2,0)));
		import(transfers.getBuffer(pos,TransferDirection(0,2,1)));
		import(transfers.getBuffer(pos,TransferDirection(0,2,2)));


		import(transfers.getBuffer(pos,TransferDirection(1,0,0)));
		import(transfers.getBuffer(pos,TransferDirection(1,0,1)));
		import(transfers.getBuffer(pos,TransferDirection(1,0,2)));

		import(transfers.getBuffer(pos,TransferDirection(1,1,0)));
		// skipped: import(transfers.getBuffer(pos,TransferDirection(1,1,1)));
		import(transfers.getBuffer(pos,TransferDirection(1,1,2)));

		import(transfers.getBuffer(pos,TransferDirection(1,2,0)));
		import(transfers.getBuffer(pos,TransferDirection(1,2,1)));
		import(transfers.getBuffer(pos,TransferDirection(1,2,2)));


		import(transfers.getBuffer(pos,TransferDirection(2,0,0)));
		import(transfers.getBuffer(pos,TransferDirection(2,0,1)));
		import(transfers.getBuffer(pos,TransferDirection(2,0,2)));

		import(transfers.getBuffer(pos,TransferDirection(2,1,0)));
		import(transfers.getBuffer(pos,TransferDirection(2,1,1)));
		import(transfers.getBuffer(pos,TransferDirection(2,1,2)));

		import(transfers.getBuffer(pos,TransferDirection(2,2,0)));
		import(transfers.getBuffer(pos,TransferDirection(2,2,1)));
		import(transfers.getBuffer(pos,TransferDirection(2,2,2)));

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
	void outputNumberOfParticlesPerCell(const Cells& cells, const std::string& outputFilename) {
		// TODO: implement output facilities for large problems
		assert_le(cells.size(), (coordinate_type{ 32,32,32 })) << "Unable to dump data for such a large cell grid at this time";


		// output particles per cell
		allscale::api::user::algorithm::async([=,&cells]() {
			auto& manager = allscale::api::core::FileIOManager::getInstance();
			auto text = manager.createEntry(outputFilename);
			auto out = manager.openOutputStream(text);

			// output dimensions
			out << cells.size() << "\n";

			std::uint64_t total = 0;
			for(std::int64_t i = 0; i < cells.size().x; ++i) {
				for(std::int64_t j = 0; j < cells.size().y; ++j) {
					for(std::int64_t k = 0; k < cells.size().z; ++k) {
						coordinate_type p{i,j,k};
						out << p.x << "," << p.y << "," << p.z << ":";
						out << cells[p].particles.size() << "\n";
						total += cells[p].particles.size();
					}
				}
			}
			out << "\nTotal: " << total << "\n";

			manager.close(out);
		}).wait();

		//allscale::api::user::algorithm::pfor(cells.size(), [&](const auto& index) {
		//	streamObject.atomic([&](auto& out) { 
		//		out << index.x << "," << index.y << "," << index.z << ":";
		//		out << cells[index].particles.size() << "\n"; });
		//});
	}

	/**
	* This function prints the positions of selected particles for visualization purposes
	*/
	template<typename StreamObject>
	void outputParticlePositions(const Cells& cells, StreamObject& out) {
		// TODO: implement output facilities for large problems
		assert_le(cells.size(), (coordinate_type{ 32,32,32 })) << "Unable to dump data for such a large cell grid at this time";

		const auto& size = cells.size();
		
		for(int i=0; i<size.x; i++) {
			for(int j=0; j<size.y; j++) {
				for(int k=0; k<size.z; k++) {
					for(const auto& p : cells[{i,j,k}].particles) {
						const auto& pos = p.position;
						out << pos.x << " " << pos.y << " " << pos.z << "\n";
					}
				}
			}
		}

	}

	/**
	* This function prints the positions of selected particles for visualization purposes
	*/
	void outputParticlePositions(const Cells& cells, const std::string& filename) {
		auto& manager = allscale::api::core::FileIOManager::getInstance();
		auto text = manager.createEntry(filename);
		auto out = manager.openOutputStream(text);
		outputParticlePositions(cells, out);
		manager.close(out);
	}

} // end namespace ipic3d
