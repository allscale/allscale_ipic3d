#pragma once

#include <vector>
#include <random>

#include "allscale/api/core/io.h"
#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/user/algorithm/async.h"
#include "allscale/api/user/algorithm/preduce.h"
#include "allscale/utils/serializer.h"
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
	 * Computes the coordinates of a cell the given particle should be located in.
	 */
	utils::Coordinate<3> getCellCoordinates(const UniverseProperties& universeProperties, const Particle& p) {
		auto cordf = elementwiseDivision((p.position - universeProperties.origin),universeProperties.cellWidth);
		return { std::int64_t(cordf.x), std::int64_t(cordf.y), std::int64_t(cordf.z) };
	}

	/**
	 * Tests whether a given particle is within the universe.
	 * @param universeProperties the properties of this universe
	 * @param p the particle to be tested
	 * @return true if the particle is inside, false otherwise
	 */
	bool isInsideUniverse(const UniverseProperties& universeProperties, const Particle& p) {
		Vector3<double> zero = universeProperties.origin;
		Vector3<double> size {
			universeProperties.origin.x + universeProperties.size.x * universeProperties.cellWidth.x,
			universeProperties.origin.y + universeProperties.size.y * universeProperties.cellWidth.y,
			universeProperties.origin.z + universeProperties.size.z * universeProperties.cellWidth.z
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
			struct electron : public allscale::utils::trivially_serializable {

				Particle operator()() const {
					Particle p;
					p.q = -1.0;
					p.qom = -25.0;
					return p;
				}

				void seed(std::uint32_t) {}

			};

			// a generator for protons
			struct proton : public allscale::utils::trivially_serializable {

				Particle operator()() const {
					Particle p;
					p.q = 1.0;
					p.qom = 1.0;
					return p;
				}

				void seed(std::uint32_t) {}
			};

		}

		namespace vector {

			// a generator for uniformly distributed vector3 instances
			class uniform {

				Vector3<double> min;
				Vector3<double> max;

				std::uniform_real_distribution<> x;
				std::uniform_real_distribution<> y;
				std::uniform_real_distribution<> z;

				std::minstd_rand randGen;

			public:

				uniform(const Vector3<double>& min, const Vector3<double>& max, std::uint32_t seed = 0)
					: min(min), max(max), x(min.x,max.x), y(min.y,max.y), z(min.z,max.z), randGen(seed) {}

				Vector3<double> operator()() {
					return {
						x(randGen), y(randGen), z(randGen)
					};
				}

				void seed(std::uint32_t seed) {
					randGen.seed(seed);
				}

				void store(allscale::utils::ArchiveWriter& out) const {
					out.write(min);
					out.write(max);
				}

				static uniform load(allscale::utils::ArchiveReader& in) {
					auto min = in.read<Vector3<double>>();
					auto max = in.read<Vector3<double>>();
					return { min, max };
				}
			};

			// a generator for uniformly distributed vector3 instances
			class uniform_r {
				Vector3<double> R1;
				Vector3<double> R2;

				std::uniform_real_distribution<> rho1;
				std::uniform_real_distribution<> rho2;
				std::uniform_real_distribution<> rho3;

				std::minstd_rand randGen;

			public:

				uniform_r(const Vector3<double>& min, const Vector3<double>& max, std::uint32_t seed = 0)
					: R1(min), R2(max), rho1(0.0,1.0), rho2(0.0,1.0), rho3(0.0,1.0), randGen(seed) {}

				Vector3<double> operator()() {
					double rh1 = rho1(randGen);
					double rh2 = rho2(randGen);
					double rh3 = rho3(randGen);
					double nu = (1.0 - 2.0 * rh2);
					return {
						pow(pow(R1.x,3)+(pow(R2.x,3)-pow(R1.x,3))*rh1,1.0/3.0) * pow(1-nu*nu, 1.0/2.0) * cos(2.0*M_PI*rh3),
						pow(pow(R1.y,3)+(pow(R2.y,3)-pow(R1.y,3))*rh1,1.0/3.0) * pow(1-nu*nu, 1.0/2.0) * sin(2.0*M_PI*rh3),
						pow(pow(R1.z,3)+(pow(R2.z,3)-pow(R1.z,3))*rh1,1.0/3.0) * nu
					};
				}

				void seed(std::uint32_t seed) {
					randGen.seed(seed);
				}

				void store(allscale::utils::ArchiveWriter& out) const {
					out.write(R1);
					out.write(R2);
				}

				static uniform_r load(allscale::utils::ArchiveReader& in) {
					auto min = in.read<Vector3<double>>();
					auto max = in.read<Vector3<double>>();
					return { min, max };
				}

			};

			// a generator for normal distributed vector3 instances
			class normal {

				Vector3<double> mean;
				Vector3<double> stddev;

				ziggurat_normal_distribution rand;

			public:

				normal(const Vector3<double>& mean, const Vector3<double>& stddev, std::uint32_t seed = 0)
					: mean(mean), stddev(stddev), rand(seed) {}

				Vector3<double> operator()() {
					return {
						mean.x + stddev.x * rand(),
						mean.y + stddev.y * rand(),
						mean.z + stddev.z * rand()
					};
				}

				void seed(std::uint32_t seed) {
					rand = ziggurat_normal_distribution(seed);
				}

				void store(allscale::utils::ArchiveWriter& out) const {
					out.write(mean);
					out.write(stddev);
				}

				static normal load(allscale::utils::ArchiveReader& in) {
					auto mean   = in.read<Vector3<double>>();
					auto stddev = in.read<Vector3<double>>();
					return { mean, stddev };
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

			void seed(std::uint32_t seed) {
				speciesGen.seed(seed);
				seed = (1023 * seed) & seed;
				posGen.seed(seed);
				seed = (1023 * seed) & seed;
				velGen.seed(seed);
			}

			void store(allscale::utils::ArchiveWriter& out) const {
				out.write(speciesGen);
				out.write(posGen);
				out.write(velGen);
			}

			static generic_particle_generator load(allscale::utils::ArchiveReader& in) {
				auto a = in.read<SpeciesDist>();
				auto b = in.read<PositionDist>();
				auto c = in.read<VelocityDist>();
				return { b, c, a };
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

			uniform(super&& base) : super(std::move(base)) {}

			void store(allscale::utils::ArchiveWriter& out) const {
				super::store(out);
			}

			static uniform load(allscale::utils::ArchiveReader& in) {
				return super::load(in);
			}
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

			normal(super&& base) : super(std::move(base)) {}

			void store(allscale::utils::ArchiveWriter& out) const {
				super::store(out);
			}

			static normal load(allscale::utils::ArchiveReader& in) {
				return super::load(in);
			}
		};

		template<typename SpeciesGen = species::electron>
		struct uniform_pos_normal_speed : public generic_particle_generator<vector::uniform,vector::normal,SpeciesGen> {

			using super = generic_particle_generator<vector::uniform,vector::normal,SpeciesGen>;

			uniform_pos_normal_speed(
					const Vector3<double>& minPos,
					const Vector3<double>& maxPos,
					const Vector3<double>& center,
					const Vector3<double>& stddev,
					std::uint32_t seed = 0
			) : super({minPos,maxPos,seed+1},{center,stddev,seed+2},SpeciesGen()) {}

			uniform_pos_normal_speed(
					const SpeciesGen& speciesGen,
					const Vector3<double>& minPos,
					const Vector3<double>& maxPos,
					const Vector3<double>& center,
					const Vector3<double>& stddev,
					std::uint32_t seed = 0
			) : super({minPos,maxPos,seed+1},{center,stddev,seed+2},speciesGen) {}

			uniform_pos_normal_speed(super&& base) : super(std::move(base)) {}

			void store(allscale::utils::ArchiveWriter& out) const {
				super::store(out);
			}

			static uniform_pos_normal_speed load(allscale::utils::ArchiveReader& in) {
				return super::load(in);
			}

		};

		template<typename SpeciesGen = species::electron>
		struct uniform_pos_normal_speed_r : public generic_particle_generator<vector::uniform_r,vector::normal,SpeciesGen> {

			using super = generic_particle_generator<vector::uniform_r,vector::normal,SpeciesGen>;

			uniform_pos_normal_speed_r(
					const Vector3<double>& minPos,
					const Vector3<double>& maxPos,
					const Vector3<double>& center,
					const Vector3<double>& stddev,
					std::uint32_t seed = 0
			) : super({minPos,maxPos,seed+1},{center,stddev,seed+2},SpeciesGen()) {}

			uniform_pos_normal_speed_r(
					const SpeciesGen& speciesGen,
					const Vector3<double>& minVel,
					const Vector3<double>& maxVel,
					const Vector3<double>& center,
					const Vector3<double>& stddev,
					std::uint32_t seed = 0
			) : super({minVel,maxVel,seed+1},{center,stddev,seed+2},speciesGen) {}

			uniform_pos_normal_speed_r(super&& base) : super(std::move(base)) {}

			void store(allscale::utils::ArchiveWriter& out) const {
				super::store(out);
			}

			static uniform_pos_normal_speed_r load(allscale::utils::ArchiveReader& in) {
				return super::load(in);
			}

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

			void seed(std::uint32_t seed) {
				dist.seed(seed);
			}

			void store(allscale::utils::ArchiveWriter& out) const {
				out.write(dist);
				out.write(center);
				out.write(radius);
			}

			static spherical load(allscale::utils::ArchiveReader& in) {
				auto dist = in.read<Distribution>();
				auto center = in.read<Vector3<double>>();
				auto radius = in.read<double>();
				return { dist, center, radius };
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
		auto gridSize = properties.size;
		Cells cells(gridSize);

		// get a private copy of the distribution generator
		auto next = dist;

		double e = 1.602176565e-19; // Elementary charge (Coulomb)
 		double m = 1.672621777e-27; // Proton mass (kg)

 		// Phase 1: approximate particle distribution

		// just some info about the progress
		std::cout << "Approximating particle distribution ...\n";

 		// create data item with distribution approximation
		auto numCells = properties.size.x * properties.size.y * properties.size.z;
		std::vector<float> distribution(numCells);
		std::vector<std::uint64_t> particleCount(numCells);

		// compute particles per cell
 		auto numPseudoParticles = numCells * 100;
		float particlesPerPseudoParticle = numParticles / float(numPseudoParticles);


		auto flatten = [=](const coordinate_type& pos) {
			return (pos.x * gridSize.y + pos.y) * gridSize.z + pos.z;
		};


 		// distribute pseudo particles
 		{

 			// compute the particle distribution approximation
 			for(int i=0; i<numPseudoParticles; i++) {
 				auto p = next();
 				while (!isInsideUniverse(properties,p)) p = next();
 				distribution[flatten(getCellCoordinates(properties,p))] += particlesPerPseudoParticle;
 			}

 			// sum up the currently distributed number of particles
 			std::uint64_t sum = 0;
 			std::uint64_t max = 0;
 			std::uint64_t min = std::numeric_limits<std::uint64_t>::max();
 			for(int i=0; i<properties.size.x; ++i) {
 				for(int j=0; j<properties.size.y; ++j) {
 					for(int k=0; k<properties.size.z; ++k) {
 						particleCount[flatten({i,j,k})] = distribution[flatten({i,j,k})];
 						sum += particleCount[flatten({i,j,k})];
 						max = std::max(max,particleCount[flatten({i,j,k})]);
 						min = std::min(min,particleCount[flatten({i,j,k})]);
 					}
 				}
 			}

 			// correct for rounding errors
 			std::int64_t missing = numParticles - sum;

 			// no error, nothing to fix
 			if (missing != 0) {

				// apply corrections
				int correct = (missing < 0) ? -((-missing / numCells) + 1) : (missing / numCells + 1);
				std::uint64_t error = std::abs(missing);
				for(int i=0; i<properties.size.x; ++i) {
					for(int j=0; j<properties.size.y; ++j) {
						for(int k=0; k<properties.size.z; ++k) {
							// correct for remaining particles (to be evenly balance)
							std::uint64_t linPos = flatten({i,j,k});
							if (linPos < error) {
								particleCount[flatten({i,j,k})] += correct;
								sum += correct;
								max = std::max(max,particleCount[flatten({i,j,k})]);
								min = std::min(min,particleCount[flatten({i,j,k})]);
							}
						}
					}
				}
 			}

 			// test that total sum of particles is correct
 			assert_eq(sum,numParticles);

 			// print seed summary
 			std::cout << "Number of particles in cells (min/avg/max): " << min << "/" << (sum/numCells) << "/" << max << "\n";
 		}


 		// Phase 2: realize approximated particle distribution

		std::cout << "Populating cells ...\n";

 		// initialize each cell in parallel
		pfor(properties.size, [=,&cells](const auto& pos) {

			// get targeted cell
			auto& cell = cells[pos];

			// get cell corners
			auto& width = properties.cellWidth;
			Vector3<double> low { width.x * pos.x, width.y * pos.y, width.z * pos.z };
			low += properties.origin;
			Vector3<double> hig = low + width;

			// compute seed for this cell
			auto seed = (uint32_t)(((pos.x * 1023) + pos.y) * 1023 + pos.z);

			// create a copy of the particle distribution and re-seed it
			auto myNext = next;
			myNext.seed(seed);

			// create a uniform position distribution for this domain
			distribution::vector::uniform next_position(low,hig,seed);

			// get number of particles to be generated in this cell
			auto localParticles = particleCount[flatten(pos)];

			// generate particles
			for(std::uint64_t i=0; i<localParticles; i++) {
				auto p = myNext();
				p.position = next_position();
				p.q = e;
				p.qom = e / m;
				cell.particles.push_back(p);
			}

			// make sure all those particles have been valid
			assert_true(verifyCorrectParticlesPositionInCell(properties,cell,pos));

		});

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
			low += properties.origin;
			Vector3<double> hig = low + width;

			// create a uniform distribution
			auto seed = (uint32_t)(((pos.x * 1023) + pos.y) * 1023 + pos.z);

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
	* This function aggregates the density contributions into the current density on the nodes
	*
	* @param universeProperties the properties of this universe
	* @param densityContributions the density contributions
	* @param pos the coordinates of this current density on the grid
	* @param density the current density output
	*/
	void aggregateDensityContributions(const UniverseProperties& universeProperties, const CurrentDensity& densityContributions, const utils::Coordinate<3>& pos, DensityNode& density) {

		auto size = densityContributions.size();
		auto curDensityContributionPos = pos * 2;

		for(int i = 0; i < 2; ++i) {
			for(int j = 0; j < 2; ++j) {
				for(int k = 0; k < 2; ++k) {
					utils::Coordinate<3> cur = curDensityContributionPos + utils::Coordinate<3>{i - 1, j - 1, k - 1};
					if(cur[0] < 0 || cur[0] >= size[0]) continue;
					if(cur[1] < 0 || cur[1] >= size[1]) continue;
					if(cur[2] < 0 || cur[2] >= size[2]) continue;

					density.J += densityContributions[cur].J;
				}
			}
		}

		const double vol = universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;
		density.J = density.J / vol / 8.0; // divide by 8 to average the value contributed by 8 neighboring cells
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

		assert_le(0,pos.x); assert_le(pos.x,1);
		assert_le(0,pos.y); assert_le(pos.y,1);
		assert_le(0,pos.z); assert_le(pos.z,1);

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
	 * @param properties the properties of this universe
	 * @param cell the cell whose particles are moved
	 * @param pos the coordinates of this cell in the grid
	 * @param field the most recently computed state of the surrounding force fields
	 */
	void moveParticles(const UniverseProperties& properties, Cell& cell, const utils::Coordinate<3>& pos, const Field& /*field*/) {

		assert_true(pos.dominatedBy(properties.size)) << "Position " << pos << " is outside universe of size " << properties.size;

		// quick-check
		if (cell.particles.empty()) return;

		// -- move the particles in space --

//		// extract forces
//		// TODO: move this to some C++ structure
//		Vector3<double> Es[2][2][2];
//		Vector3<double> Bs[2][2][2];
//		for(int i=0; i<2; i++) {
//			for(int j=0; j<2; j++) {
//				for(int k=0; k<2; k++) {
//					utils::Coordinate<3> cur({pos[0]+i+1,pos[1]+j+1,pos[2]+k+1});
//					Es[i][j][k] = field[cur].E;
//					Bs[i][j][k] = field[cur].B;
//				}
//			}
//		}

		//const auto cellOrigin = getOriginOfCell(pos, properties);

		//double vol = properties.cellWidth.x * properties.cellWidth.y * properties.cellWidth.z;

		double magneticFieldTemp = -properties.externalMagneticField.z * pow(properties.planetRadius, 3);

		// update particles
//		allscale::api::user::algorithm::pfor(cell.particles, [&](Particle& p){
		for(std::size_t i = 0; i < cell.particles.size(); ++i) {
			Particle& p = cell.particles[i];
			// Docu: https://www.particleincell.com/2011/vxb-rotation/
			// Code: https://www.particleincell.com/wp-content/uploads/2011/07/ParticleIntegrator.java

//			// get the fractional distance of the particle from the cell origin
//			const auto relPos = allscale::utils::elementwiseDivision((p.position - cellOrigin), (properties.cellWidth));
//
//			// interpolate
//			auto E = trilinearInterpolationF2P(Es, relPos, vol);
//			auto B = trilinearInterpolationF2P(Bs, relPos, vol);

			// calculate 3 Cartesian components of the magnetic field
			double fac1 =  magneticFieldTemp / pow(allscale::utils::sumOfSquares(p.position), 2.5);
			Vector3<double> E, B;
			E = {0.0, 0.0, 0.0};
			B.x = 3.0 * p.position.x * p.position.z * fac1;
			B.y = 3.0 * p.position.y * p.position.z * fac1;
			B.z = (2.0 * pow(p.position.z, 2) - pow(p.position.x, 2) - pow(p.position.y, 2)) * fac1;

			// adaptive sub-cycling for computing velocity
			double B_mag = allscale::utils::sumOfSquares(B);
			double dt_sub = M_PI * properties.speedOfLight / (4.0 * fabs(p.qom) * B_mag);
			int sub_cycles = int(properties.dt / dt_sub) + 1;
			sub_cycles = std::min(sub_cycles, 100);
			dt_sub = properties.dt / double(sub_cycles);

			for (int cyc_cnt = 0; cyc_cnt < sub_cycles; cyc_cnt++) {
				// update velocity
				p.updateVelocity(E, B, dt_sub);

				// update position
				p.updatePosition(dt_sub);
			}
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

			neighbors[0][0][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, -1, -1 } +size) % size, TransferDirection{ 2, 2, 2 });
			neighbors[0][0][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, -1, 0 } +size) % size, TransferDirection{ 2, 2, 1 });
			neighbors[0][0][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, -1, 1 } +size) % size, TransferDirection{ 2, 2, 0 });

			neighbors[0][1][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, 0, -1 } +size) % size, TransferDirection{ 2, 1, 2 });
			neighbors[0][1][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, 0, 0 } +size) % size, TransferDirection{ 2, 1, 1 });
			neighbors[0][1][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, 0, 1 } +size) % size, TransferDirection{ 2, 1, 0 });

			neighbors[0][2][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, 1, -1 } +size) % size, TransferDirection{ 2, 0, 2 });
			neighbors[0][2][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, 1, 0 } +size) % size, TransferDirection{ 2, 0, 1 });
			neighbors[0][2][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{ -1, 1, 1 } +size) % size, TransferDirection{ 2, 0, 0 });


			neighbors[1][0][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, -1, -1 } +size) % size, TransferDirection{ 1, 2, 2 });
			neighbors[1][0][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, -1, 0 } +size) % size, TransferDirection{ 1, 2, 1 });
			neighbors[1][0][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, -1, 1 } +size) % size, TransferDirection{ 1, 2, 0 });

			neighbors[1][1][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, 0, -1 } +size) % size, TransferDirection{ 1, 1, 2 });
			neighbors[1][1][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, 0, 0 } +size) % size, TransferDirection{ 1, 1, 1 });
			neighbors[1][1][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, 0, 1 } +size) % size, TransferDirection{ 1, 1, 0 });

			neighbors[1][2][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, 1, -1 } +size) % size, TransferDirection{ 1, 0, 2 });
			neighbors[1][2][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, 1, 0 } +size) % size, TransferDirection{ 1, 0, 1 });
			neighbors[1][2][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  0, 1, 1 } +size) % size, TransferDirection{ 1, 0, 0 });


			neighbors[2][0][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, -1, -1 } +size) % size, TransferDirection{ 0, 2, 2 });
			neighbors[2][0][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, -1, 0 } +size) % size, TransferDirection{ 0, 2, 1 });
			neighbors[2][0][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, -1, 1 } +size) % size, TransferDirection{ 0, 2, 0 });

			neighbors[2][1][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, 0, -1 } +size) % size, TransferDirection{ 0, 1, 2 });
			neighbors[2][1][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, 0, 0 } +size) % size, TransferDirection{ 0, 1, 1 });
			neighbors[2][1][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, 0, 1 } +size) % size, TransferDirection{ 0, 1, 0 });

			neighbors[2][2][0] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, 1, -1 } +size) % size, TransferDirection{ 0, 0, 2 });
			neighbors[2][2][1] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, 1, 0 } +size) % size, TransferDirection{ 0, 0, 1 });
			neighbors[2][2][2] = &transfers.getBuffer((pos + utils::Coordinate<3>{  1, 1, 1 } +size) % size, TransferDirection{ 0, 0, 0 });

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
			//			allscale::api::user::algorithm::pfor(std::size_t(0),cell.particles.size(),[&](std::size_t index){
			for(std::size_t index = 0; index<cell.particles.size(); ++index) {

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

				// remove particles from inside the sphere
				auto diff = p.position - universeProperties.objectCenter;
				double r2 = allscale::utils::sumOfSquares(diff);
				if(r2 <= universeProperties.planetRadius * universeProperties.planetRadius) {
					continue;
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
				if(targets[i]) {
					targets[i]->push_back(cell.particles[i]);
				}
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

        std::array<std::vector<Particle>*, 26> buffers;

        buffers[ 0] = &transfers.getBuffer(pos,TransferDirection(0,0,0));
        buffers[ 1] = &transfers.getBuffer(pos,TransferDirection(0,0,1));
        buffers[ 2] = &transfers.getBuffer(pos,TransferDirection(0,0,2));

        buffers[ 3] = &transfers.getBuffer(pos,TransferDirection(0,1,0));
        buffers[ 4] = &transfers.getBuffer(pos,TransferDirection(0,1,1));
        buffers[ 5] = &transfers.getBuffer(pos,TransferDirection(0,1,2));

        buffers[ 6] = &transfers.getBuffer(pos,TransferDirection(0,2,0));
        buffers[ 7] = &transfers.getBuffer(pos,TransferDirection(0,2,1));
        buffers[ 8] = &transfers.getBuffer(pos,TransferDirection(0,2,2));

        buffers[ 9] = &transfers.getBuffer(pos,TransferDirection(1,0,0));
        buffers[10] = &transfers.getBuffer(pos,TransferDirection(1,0,1));
        buffers[11] = &transfers.getBuffer(pos,TransferDirection(1,0,2));

        buffers[12] = &transfers.getBuffer(pos,TransferDirection(1,1,0));
		// skipped: import(transfers.getBuffer(pos,TransferDirection(1,1,1)));
        buffers[13] = &transfers.getBuffer(pos,TransferDirection(1,1,2));

        buffers[14] = &transfers.getBuffer(pos,TransferDirection(1,2,0));
        buffers[15] = &transfers.getBuffer(pos,TransferDirection(1,2,1));
        buffers[16] = &transfers.getBuffer(pos,TransferDirection(1,2,2));

        buffers[17] = &transfers.getBuffer(pos,TransferDirection(2,0,0));
        buffers[18] = &transfers.getBuffer(pos,TransferDirection(2,0,1));
        buffers[19] = &transfers.getBuffer(pos,TransferDirection(2,0,2));

        buffers[20] = &transfers.getBuffer(pos,TransferDirection(2,1,0));
        buffers[21] = &transfers.getBuffer(pos,TransferDirection(2,1,1));
        buffers[22] = &transfers.getBuffer(pos,TransferDirection(2,1,2));

        buffers[23] = &transfers.getBuffer(pos,TransferDirection(2,2,0));
        buffers[24] = &transfers.getBuffer(pos,TransferDirection(2,2,1));
        buffers[25] = &transfers.getBuffer(pos,TransferDirection(2,2,2));

        std::size_t newSize = 0;
        for (auto buffer: buffers)
            newSize += buffer->size();
        cell.particles.reserve(cell.particles.size() + newSize);

		// along all 26 directions (center is not relevant)
		import(*buffers[ 0]);
		import(*buffers[ 1]);
		import(*buffers[ 2]);

		import(*buffers[ 3]);
		import(*buffers[ 4]);
		import(*buffers[ 5]);

		import(*buffers[ 6]);
		import(*buffers[ 7]);
		import(*buffers[ 8]);

		import(*buffers[ 9]);
		import(*buffers[10]);
		import(*buffers[11]);

		import(*buffers[12]);
		// skipped: import(transfers.getBuffer(pos,TransferDirection(1,1,1)));
		import(*buffers[13]);

		import(*buffers[14]);
		import(*buffers[15]);
		import(*buffers[16]);

		import(*buffers[17]);
		import(*buffers[18]);
		import(*buffers[19]);

		import(*buffers[20]);
		import(*buffers[21]);
		import(*buffers[22]);

		import(*buffers[23]);
		import(*buffers[24]);
		import(*buffers[25]);

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
