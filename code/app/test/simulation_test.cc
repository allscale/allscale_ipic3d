#include <gtest/gtest.h>

#include "ipic3d/app/simulator.h"
#include "ipic3d/app/universe.h"
#include "ipic3d/app/common.h"

#include "allscale/api/user/operator/pfor.h"

namespace ipic3d {


	namespace detail {


		template<typename Mover>
		void testSingleParticle() {

			// this test checks whether particles are properly migrated between cells

			// Set universe properties
			UniverseProperties properties;
			properties.size = { 2,1,1 };
			properties.cellWidth = { 1,1,1 };
			properties.dt = 1;
			properties.useCase = UseCase::Test;

			// Create a universe with these properties
			Universe universe = Universe(properties);

			Field& field = universe.field;
			decltype(field.size()) zero = 0;
			allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
				field[pos].E = { 0, 0, 0 };
				field[pos].B = { 0, 0, 0 };
			});

			Cell& a = universe.cells[{0,0,0}];
			Cell& b = universe.cells[{1,0,0}];

			// insert one particle
			Particle p;
			p.position.x = p.position.y = p.position.z = 0.5;
			p.velocity.x = p.velocity.y = p.velocity.z = 0.0;
			p.mass = 1;
			p.q = 1;

			// add test particle to first cell
			a.particles.push_back(p);

			// check particle position
			ASSERT_FALSE(a.particles.empty());
			ASSERT_TRUE(b.particles.empty());

			EXPECT_EQ(0.5, a.particles.front().position.x);
			EXPECT_EQ(0.5, a.particles.front().position.y);
			EXPECT_EQ(0.5, a.particles.front().position.z);

			EXPECT_EQ(0.0, a.particles.front().velocity.x);
			EXPECT_EQ(0.0, a.particles.front().velocity.y);
			EXPECT_EQ(0.0, a.particles.front().velocity.z);

			EXPECT_EQ(1.0, a.particles.front().mass);
			EXPECT_EQ(1.0, a.particles.front().q);

			// run one simulation step
			simulateStep<detail::default_particle_to_field_projector,detail::default_field_solver,Mover>(universe);

			// check particle position
			ASSERT_FALSE(a.particles.empty());
			ASSERT_TRUE(b.particles.empty());

			EXPECT_EQ(0.5, a.particles.front().position.x);
			EXPECT_EQ(0.5, a.particles.front().position.y);
			EXPECT_EQ(0.5, a.particles.front().position.z);

			EXPECT_EQ(0.0, a.particles.front().velocity.x);
			EXPECT_EQ(0.0, a.particles.front().velocity.y);
			EXPECT_EQ(0.0, a.particles.front().velocity.z);

			EXPECT_EQ(1.0, a.particles.front().mass);
			EXPECT_EQ(1.0, a.particles.front().q);

			// change velocity and send in x direction
			Particle& p2 = a.particles.front();
			p2.velocity.x = 1.0;
			p2.velocity.y = p2.velocity.z = 0.0;

			EXPECT_EQ(0.5, a.particles.front().position.x);
			EXPECT_EQ(0.5, a.particles.front().position.y);
			EXPECT_EQ(0.5, a.particles.front().position.z);

			EXPECT_EQ(1.0, a.particles.front().velocity.x);
			EXPECT_EQ(0.0, a.particles.front().velocity.y);
			EXPECT_EQ(0.0, a.particles.front().velocity.z);

			EXPECT_EQ(1.0, a.particles.front().mass);
			EXPECT_EQ(1.0, a.particles.front().q);

			// run one simulation step
			simulateStep<detail::default_particle_to_field_projector, detail::default_field_solver, Mover>(universe);

			ASSERT_TRUE(a.particles.empty());
			ASSERT_FALSE(b.particles.empty());

			// -- check the particle position
			EXPECT_EQ(1.5,b.particles.front().position.x);
			EXPECT_EQ(0.5,b.particles.front().position.y);
			EXPECT_EQ(0.5,b.particles.front().position.z);

			EXPECT_EQ(1.0,b.particles.front().velocity.x);
			EXPECT_EQ(0.0,b.particles.front().velocity.y);
			EXPECT_EQ(0.0,b.particles.front().velocity.z);

			// change velocity and send back
			b.particles.front().velocity.x = -1;

			int niter = 1;

			simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,Mover>(niter,universe);

			ASSERT_FALSE(a.particles.empty());
			ASSERT_TRUE(b.particles.empty());

			EXPECT_EQ(0.5,a.particles.front().position.x);
			EXPECT_EQ(0.5,a.particles.front().position.y);
			EXPECT_EQ(0.5,a.particles.front().position.z);

			EXPECT_EQ(-1.0,a.particles.front().velocity.x);
			EXPECT_EQ( 0.0,a.particles.front().velocity.y);
			EXPECT_EQ( 0.0,a.particles.front().velocity.z);

		}

	}


	TEST(Simulation, SingleParticleFirstOrder) {

		// this test checks whether particles are properly migrated between cells
		detail::testSingleParticle<detail::default_particle_mover>();

	}


	TEST(Simulation, SingleParticleBorisMover) {

		// Set universe properties
		UniverseProperties properties;
		properties.size = {1,1,1};
		properties.cellWidth = { 1,1,1 };
		properties.dt = 0.1;
		properties.useCase = UseCase::Test;

		// number of steps
		unsigned niter = 9;

		// Create a universe with these properties
		Universe universe = Universe(properties);

		// configure the cell
		Cell& cell = universe.cells[{0,0,0}];

		// initialize the field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.2, 0.0, 0.0 };
			field[pos].B = { 0.2, 0.0, 0.0 };
		});

		// add one particle
		Particle p;
		p.position.x = p.position.y = 0.5;
		p.position.z = 0;

		p.velocity.x = p.velocity.y = 0;
		p.velocity.z = 1;

		p.q = 1;
		p.mass = 1;

		cell.particles.push_back(p);

		// run the simulation
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check where particle ended up
		ASSERT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected
		EXPECT_NEAR( res.position.x, 0.590, 0.001);
		EXPECT_NEAR( res.position.y, 0.589, 0.001);
		EXPECT_NEAR( res.position.z, 0.894, 0.001);

	}


	TEST(Simulation, SingleParticleBorisMoverLarmorRadius) {

		// Set universe properties
		UniverseProperties properties;
		properties.useCase = UseCase::Test;
		properties.size = { 1,1,1 };
		properties.cellWidth = { 1e4,1e4,1e4 };
		properties.dt = 3e-11;
		properties.origin = { -10.0, -10.0, -10.0 };
		properties.FieldOutputCycle = 1e6;

		// Create Universe with these properties
		Universe universe = Universe(properties);

		Cell& cell = universe.cells[{0, 0, 0}];

		// initialize field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.0, 0.0, 0.0  };
			field[pos].B = { 0.0, 0.0, 0.01 };
		});

		utils::Coordinate<3> field_pos = { 1, 1, 1 };

		// add one particle
		Particle p;
		p.position.x = p.position.y = p.position.z = 0.0;

		p.velocity.x = p.velocity.z = 0.0;
		p.velocity.y = 1e5;

		p.q = -1.602e-19;
		p.mass = 9.109e-31;

		// compute Larmor radius
		double rL = p.mass * p.velocity.y / (fabs(p.q) * field[field_pos].B.z);
		EXPECT_NEAR( rL, 5.686e-05, 1e-06 );

		// re-initialize the x coordinate
		p.position.x = rL;

		// push velocity back in time by 1/2 dt
		// 		this is purely done to compare against the Matlab version 
		p.updateVelocityBorisStyle(field[field_pos].E, field[field_pos].B, -0.5*properties.dt);

		// number of steps
		unsigned numSteps = 1000;
		// run the simulation
		//simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(niter,universe);
		for(unsigned i = 0; i < numSteps; ++i) {
			auto E = field[field_pos].E;
			auto B = field[field_pos].B;

			// update velocity
			p.updateVelocityBorisStyle(E, B, properties.dt);

			// update position
			p.updatePosition(properties.dt);
		}
		

		// check where particle ended up
		cell.particles.push_back(p);
		ASSERT_FALSE(cell.particles.empty());

		// check that the position is close to what is expected
		// comparing against the matlab code after 10 iterations
		EXPECT_NEAR( p.position.x, -4.50134e-05, 1e-10);
		EXPECT_NEAR( p.position.y, 3.47983e-05, 1e-09);
		EXPECT_NEAR( p.position.z, 0.0,		 1e-10);

		// check that the velocity is close to what is expected
		EXPECT_NEAR( p.velocity.x, -63242.7, 1e1);
		EXPECT_NEAR( p.velocity.y, -77462.0, 1e1);
		EXPECT_NEAR( p.velocity.z, 0.0,     1e-02);
	}


	TEST(Simulation, SingleParticleLarmorRadius) {

		// Set universe properties
		UniverseProperties properties;
		properties.useCase = UseCase::Test;
		properties.size = { 1,1,1 };
		properties.cellWidth = { 1e4,1e4,1e4 };
		properties.dt = 3e-11;
		properties.origin = { -10.0, -10.0, -10.0 };
		properties.FieldOutputCycle = 1e6;

		// Create Universe with these properties
		Universe universe = Universe(properties);

		Cell& cell = universe.cells[{0, 0, 0}];

		// initialize field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.0, 0.0, 0.0  };
			field[pos].B = { 0.0, 0.0, 0.01 };
		});

		utils::Coordinate<3> field_pos = { 1, 1, 1 };

		// add one particle
		Particle p;
		p.position.x = p.position.y = p.position.z = 0.0;

		p.velocity.x = p.velocity.z = 0.0;
		p.velocity.y = 1e5;

		p.q = -1.602e-19;
		p.mass = 9.109e-31;

		// compute Larmor radius
		double rL = p.mass * p.velocity.y / (fabs(p.q) * field[field_pos].B.z);
		EXPECT_NEAR( rL, 5.686e-05, 1e-06 );

		// re-initialize the x coordinate
		p.position.x = rL;

		// push velocity back in time by 1/2 dt
		// 		this is purely done to compare against the Matlab version 
		p.updateVelocityBorisStyle(field[field_pos].E, field[field_pos].B, -0.5*properties.dt);

		cell.particles.push_back(p);

		// check particle position
		EXPECT_FALSE(cell.particles.empty());

		// number of steps
		unsigned numSteps = 10;
		// run the simulation
		simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(numSteps,universe);

		// check where particle ended up
		ASSERT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected
		// comparing against the matlab code after 10 iterations
		EXPECT_NEAR( res.position.x, 5.7651e-05, 1e-09);
	    EXPECT_NEAR( res.position.y, 2.9990e-05, 1e-09); 
        EXPECT_NEAR( res.position.z, 0.0, 1e-10);

		// check that the velocity is close to what is expected
		EXPECT_NEAR( p.velocity.x, 2637.59, 1e1);
		EXPECT_NEAR( p.velocity.y, 99965.2, 1e1);
		EXPECT_NEAR( p.velocity.z, 0.0,     1e-02);
	}


	TEST(Simulation, SingleParticleBorisMoverExBdrift) {

		// Set universe properties
		UniverseProperties properties;
		properties.size = { 1,1,1 };
		properties.cellWidth = { 1e4,1e4,1e4 };
		properties.dt = 0.01;
		properties.useCase = UseCase::Test;
		properties.FieldOutputCycle = 1e6;

		// Create Universe with these properties
		Universe universe = Universe(properties);

		Cell& cell = universe.cells[{0, 0, 0}];

		// initialize field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.2, 0.0, 0.0 };
			field[pos].B = { 0.0, 0.0, 1.0 };
		});

		// add one particle
		Particle p;
		p.position.x = p.position.y = p.position.z = 0.0;

		p.velocity.y = p.velocity.z = 0.0;
		p.velocity.x = 0.5;

		p.q = 1.0;
		p.mass = 1.0;

		// push velocity back in time by 1/2 dt
		// 		this is purely done to compare against the Matlab version 
		//p.updateVelocityBorisStyle(field[{0,0,0}].E, field[{0,0,0}].B, -0.5*properties.dt);

		cell.particles.push_back(p);

		// check particle position
		EXPECT_FALSE(cell.particles.empty());

		// run the simulation
		// number of steps
		unsigned niter = 5000;
		simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(niter,universe);
//		utils::Coordinate<3> pos = { 0, 0, 0 };
//		std::cout << cell.particles.front();
//		moveParticlesBorisStyle(universe.properties, cell, pos, field);
//		std::cout << cell.particles.front();
//		moveParticlesBorisStyle(universe.properties, cell, pos, field);
//		std::cout << cell.particles.front();
//		moveParticlesBorisStyle(universe.properties, cell, pos, field);
//		std::cout << cell.particles.front();
//		moveParticlesBorisStyle(universe.properties, cell, pos, field);
//		std::cout << cell.particles.front();
//		moveParticlesBorisStyle(universe.properties, cell, pos, field);
//		std::cout << cell.particles.front();
		return;

		// check where particle ended up
		ASSERT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected
		// comparing against the matlab code after 10 iterations
		EXPECT_NEAR( res.position.x, 0.005009, 1e-4 );
		EXPECT_NEAR( res.position.y, 1.0e-04 * -0.2503, 1e-06 );
		EXPECT_NEAR( res.position.z, 0.0,		 1e-10 );

		// check that the velocity is close to what is expected
		EXPECT_NEAR( res.velocity.x, 0.50197, 1e-04);
		EXPECT_NEAR( res.velocity.y, -0.005009, 1e-04);
		EXPECT_NEAR( res.velocity.z, 0.0,     1e-02);
	}


	TEST(SimulationTest, ParticleMigrationOneDirectionOneParticle) {

		// this test checks whether particles are properly migrated between cells
		// especially for the periodic boundary conditions when particles exit the domain

		// Set universe properties
		UniverseProperties properties;
		properties.size = {3,3,3};
		properties.cellWidth = { .5,.5,.5 };
		properties.dt = 0.1;
		properties.useCase = UseCase::Test;

		// Create a universe with these properties
		Universe universe = Universe(properties);

		// configure the cell
		Cell& a = universe.cells[{0,1,0}];
		Cell& b = universe.cells[{1,1,0}];
		Cell& c = universe.cells[{2,1,0}];

		// initialize the field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.2, 0.0, 0.0 };
			field[pos].B = { 0.2, 0.0, 0.0 };
		});

		// add one particle
		Particle p;
		p.position.z = 0.0; p.position.y = 0.8;
		p.position.x = 0.4;
		p.velocity.z = p.velocity.y = 0.0;
		p.velocity.x = 1.0;
		p.q = 1.0;
		p.mass = 1.0;

		// add test particle to first cell
		a.particles.push_back(p);

		// check particle position
		ASSERT_FALSE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());
		ASSERT_TRUE(c.particles.empty());


		// run the simulation
		// number of steps
		unsigned niter = 1;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check whether the particle was moved from one cell to another
		ASSERT_TRUE(a.particles.empty());
		ASSERT_FALSE(b.particles.empty());
		ASSERT_TRUE(c.particles.empty());

		EXPECT_NEAR(b.particles.front().position.x, 0.516, 1e-3);
		EXPECT_NEAR(b.particles.front().position.y, 0.8, 1e-1);
		EXPECT_NEAR(b.particles.front().position.z, 0.0, 1e-0);

		EXPECT_NEAR(b.particles.front().velocity.x, 1.16, 1e-2);
		EXPECT_NEAR(b.particles.front().velocity.y, 0.0, 1e-0);
		EXPECT_NEAR(b.particles.front().velocity.z, 0.0, 1e-0);

		// check number of particles in the domain
		int total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(1, total_particles);
	   
		// verify the proper placement of particles 
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {
			ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, universe.cells[pos], pos) );
		});


		// run the simulation to propagate the particle further
		niter = 4;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check whether the particle was moved out of the domain
		ASSERT_TRUE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());
		ASSERT_FALSE(c.particles.empty());

		EXPECT_NEAR(c.particles.front().position.x, 1.14, 1e-2);
		EXPECT_NEAR(c.particles.front().position.y, 0.8, 1e-1);
		EXPECT_NEAR(c.particles.front().position.z, 0.0, 1e-0);

		EXPECT_NEAR(c.particles.front().velocity.x, 1.8, 1e-1);
		EXPECT_NEAR(c.particles.front().velocity.y, 0.0, 1e-0);
		EXPECT_NEAR(c.particles.front().velocity.z, 0.0, 1e-0);

		// check number of particles in the domain
		total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(1, total_particles);
	   
		// verify the proper placement of particles 
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {
			ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, universe.cells[pos], pos) );
		});
	

		// run the simulation to push the particle outside the domain
		niter = 2;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check whether the particle was moved out of the domain
		ASSERT_FALSE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());
		ASSERT_TRUE(c.particles.empty());

		EXPECT_NEAR(a.particles.front().position.x, 0.048, 1e-3);
		EXPECT_NEAR(a.particles.front().position.y, 0.8, 1e-1);
		EXPECT_NEAR(a.particles.front().position.z, 0.0, 1e-0);

		EXPECT_NEAR(a.particles.front().velocity.x, 2.12, 1e-2);
		EXPECT_NEAR(a.particles.front().velocity.y, 0.0, 1e-0);
		EXPECT_NEAR(a.particles.front().velocity.z, 0.0, 1e-0);

		// check number of particles in the domain
		total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(1, total_particles);
	   
		// verify the proper placement of particles 
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {
			ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, universe.cells[pos], pos) );
		});


		// now we push the same particle in the oposite direction 
		Particle& p2 = a.particles.front();
		p2.velocity.x  = -1.0;
		niter = 1;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		EXPECT_NEAR(c.particles.front().position.x, 1.464, 1e-3);
		EXPECT_NEAR(c.particles.front().position.y, 0.8, 1e-1);
		EXPECT_NEAR(c.particles.front().position.z, 0.0, 1e-0);

		EXPECT_NEAR(c.particles.front().velocity.x, -0.84, 1e-2);
		EXPECT_NEAR(c.particles.front().velocity.y, 0.0, 1e-0);
		EXPECT_NEAR(c.particles.front().velocity.z, 0.0, 1e-0);

		// check number of particles in the domain
		total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(1, total_particles);

		// check whether the particle was moved out of the domain
		ASSERT_TRUE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());
		ASSERT_FALSE(c.particles.empty());
	}


	TEST(SimulationTest, ParticleMigrationTwoDirectionsOneParticle) {

		// this test checks whether particles are properly migrated between cells
		// especially for the periodic boundary conditions when particles exit the domain

		// Set universe properties
		UniverseProperties properties;
		properties.size = {3,3,3};
		properties.cellWidth = { .5,.5,.5 };
		properties.dt = 0.1;
		properties.useCase = UseCase::Test;

		// Create a universe with these properties
		Universe universe = Universe(properties);

		// configure the cell
		Cell& a = universe.cells[{0,0,1}];
		Cell& b = universe.cells[{1,1,1}];
		Cell& c = universe.cells[{2,2,1}];

		// initialize the field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.2, 0.2, 0.0 };
			field[pos].B = { 0.2, 0.2, 0.0 };
		});

		// add one particle
		Particle p;
		p.position.x = p.position.y = 0.4;
		p.position.z = 0.8;
		p.velocity.x = p.velocity.y = 1.0;
		p.velocity.z = 0.0;
		p.q = 1.0;
		p.mass = 1.0;

		// add test particle to first cell
		a.particles.push_back(p);

		// check particle position
		ASSERT_FALSE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());
		ASSERT_TRUE(c.particles.empty());


		// run the simulation
		// number of steps
		unsigned niter = 1;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check whether the particle was moved from one cell to another
		ASSERT_TRUE(a.particles.empty());
		ASSERT_FALSE(b.particles.empty());
		ASSERT_TRUE(c.particles.empty());

		EXPECT_NEAR(b.particles.front().position.x, 0.516, 1e-3);
		EXPECT_NEAR(b.particles.front().position.y, 0.516, 1e-3);
		EXPECT_NEAR(b.particles.front().position.z, 0.8, 1e-1);

		EXPECT_NEAR(b.particles.front().velocity.x, 1.16, 1e-2);
		EXPECT_NEAR(b.particles.front().velocity.y, 1.16, 1e-2);
		EXPECT_NEAR(b.particles.front().velocity.z, 0.0, 1e-0);

		// check number of particles in the domain
		int total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(1, total_particles);
	   
		// verify the proper placement of particles 
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {
			ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, universe.cells[pos], pos) );
		});


		// run the simulation to propagate the particle further
		niter = 4;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check whether the particle was moved out of the domain
		ASSERT_TRUE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());
		ASSERT_FALSE(c.particles.empty());

		EXPECT_NEAR(c.particles.front().position.x, 1.14, 1e-2);
		EXPECT_NEAR(c.particles.front().position.y, 1.14, 1e-2);
		EXPECT_NEAR(c.particles.front().position.z, 0.8, 1e-1);

		EXPECT_NEAR(c.particles.front().velocity.x, 1.8, 1e-1);
		EXPECT_NEAR(c.particles.front().velocity.y, 1.8, 1e-1);
		EXPECT_NEAR(c.particles.front().velocity.z, 0.0, 1e-0);

		// check number of particles in the domain
		total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(1, total_particles);
	   
		// verify the proper placement of particles 
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {
			ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, universe.cells[pos], pos) );
		});
	

		// run the simulation to push the particle outside the domain
		niter = 2;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check whether the particle was moved out of the domain
		ASSERT_FALSE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());
		ASSERT_TRUE(c.particles.empty());

		EXPECT_NEAR(a.particles.front().position.x, 0.048, 1e-3);
		EXPECT_NEAR(a.particles.front().position.y, 0.048, 1e-3);
		EXPECT_NEAR(a.particles.front().position.z, 0.8, 1e-1);

		EXPECT_NEAR(a.particles.front().velocity.x, 2.12, 1e-2);
		EXPECT_NEAR(a.particles.front().velocity.y, 2.12, 1e-2);
		EXPECT_NEAR(a.particles.front().velocity.z, 0.0, 1e-0);

		// check number of particles in the domain
		total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(1, total_particles);
	   
		// verify the proper placement of particles 
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {
			ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, universe.cells[pos], pos) );
		});


		// now we push the same particle in the oposite direction 
		Particle& p2 = a.particles.front();
		p2.velocity.x = p2.velocity.y = -1.0;
		niter = 1;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		EXPECT_NEAR(c.particles.front().position.x, 1.464, 1e-3);
		EXPECT_NEAR(c.particles.front().position.y, 1.464, 1e-3);
		EXPECT_NEAR(c.particles.front().position.z, 0.8, 1e-1);

		EXPECT_NEAR(c.particles.front().velocity.x, -0.84, 1e-2);
		EXPECT_NEAR(c.particles.front().velocity.y, -0.84, 1e-2);
		EXPECT_NEAR(c.particles.front().velocity.z, 0.0, 1e-0);

		// check number of particles in the domain
		total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(1, total_particles);

		// check whether the particle was moved out of the domain
		ASSERT_TRUE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());
		ASSERT_FALSE(c.particles.empty());
	}


	TEST(SimulationTest, ParticleMigrationManyParticles) {

		// this test checks whether particles are properly migrated between cells

		// Set universe properties
		UniverseProperties properties;
		properties.size = {4,4,4};
		properties.cellWidth = { .5,.5,.5 };
		properties.dt = 0.1;
		properties.useCase = UseCase::Test;

		// Create a universe with these properties
		Universe universe = Universe(properties);

		// configure the cell
		Cell& a = universe.cells[{0,0,0}];
		Cell& b = universe.cells[{1,0,0}];
		Cell& c = universe.cells[{0,1,0}];
		Cell& d = universe.cells[{1,1,0}];
		Cell& e = universe.cells[{0,0,1}];
		Cell& f = universe.cells[{1,0,1}];
		Cell& g = universe.cells[{0,1,1}];
		Cell& h = universe.cells[{1,1,1}];

		// initialize the field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.2, 0.2, 0.2 };
			field[pos].B = { 0.2, 0.2, 0.2 };
		});

		// add one particle
		Particle p0, p1, p2, p3, p4, p5, p6, p7;
		p0.position.x = 0.4; p0.position.y = p0.position.z = 0.4;
		p0.velocity.x = 0.6; p0.velocity.y = p0.velocity.z = 0.6;
		p0.q = p0.mass = 1;
		p1.position.x = 0.8; p1.position.y = p1.position.z = 0.4;
		p1.velocity.x = 0.3; p1.velocity.y = p1.velocity.z = 0.3;
		p1.q = p1.mass = 1;
		p2.position.y = 0.8; p2.position.x = p2.position.z = 0.4;
		p2.velocity.x = 0.1; p2.velocity.y = p2.velocity.z = 0.1;
		p2.q = p2.mass = 1;
		p3.position.x = 0.75; p3.position.y = 0.75; p3.position.z = 0.4;
		p3.velocity.x = 0.8; p3.velocity.y = p3.velocity.z = 0.8;
		p3.q = p3.mass = 1;
		p4 = p0;
		p4.position.z = 0.6;
		p5 = p1;
		p5.position.z = 0.8;
		p6 = p3;
		p6.position.x = 0.25;
		p6.position.z = 0.75;
		p7 = p3;
		p7.position.z = 0.8;

		// add test particle to first cell
		a.particles.push_back(p0);
		b.particles.push_back(p1);
		c.particles.push_back(p2);
		d.particles.push_back(p3);
		e.particles.push_back(p4);
		f.particles.push_back(p5);
		g.particles.push_back(p6);
		h.particles.push_back(p7);

		// check number of particles in the domain
		int total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(8, total_particles);
	   
		// verify the proper placement of particles 
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {
			ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, universe.cells[pos], pos) );
		});

		// run the simulation
		// number of steps
		unsigned niter = 20;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check number of particles in the domain
		total_particles = countParticlesInDomain(universe);
		EXPECT_EQ(8, total_particles);
	   
		// verify the proper placement of particles 
		allscale::api::user::pfor(zero, properties.size, [&](const utils::Coordinate<3>& pos) {
			ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, universe.cells[pos], pos) );
		});
	}

} // end namespace ipic3d
