#include <gtest/gtest.h>

#include "ipic3d/app/simulator.h"
#include "ipic3d/app/universe.h"

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
		properties.size = { 1,1,1 };
		properties.cellWidth = { 1e4,1e4,1e4 };
		properties.dt = 3e-11;
		properties.useCase = UseCase::Test;

		// number of steps
		unsigned niter = 10;

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

		// add one particle
		Particle p;
		p.position.x = p.position.y = p.position.z = 0.0;

		p.velocity.x = p.velocity.z = 0.0;
		p.velocity.y = 1e5;

		p.q = -1.602e-19;
		p.mass = 9.109e-31;

		// compute Larmor radius
		double rL = p.mass * p.velocity.y / (fabs(p.q) * field[{0,0,0}].B.z);
		EXPECT_NEAR( rL, 5.686e-05, 1e-06 );
		p.position.x = rL;

		// push velocity back in time by 1/2 dt
		p.updateVelocityBorisStyle(field[{0,0,0}].E, field[{0,0,0}].B, -0.5*properties.dt);

		cell.particles.push_back(p);

		// check particle position
		EXPECT_FALSE(cell.particles.empty());

		// run the simulation
		simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(niter,universe);

		// check where particle ended up
		ASSERT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected
		// comparing against the matlab code after 10 iterations
		EXPECT_NEAR( res.position.x, 5.7652e-05, 1e-06);
		EXPECT_NEAR( res.position.y, 2.9990e-05, 1e-06);
		EXPECT_NEAR( res.position.z, 0.0,		 1e-06);

		// check that the velocity is close to what is expected
		EXPECT_NEAR( res.velocity.x, 2637.59, 1e2);
		EXPECT_NEAR( res.velocity.y, 99965.2, 1e2);
		EXPECT_NEAR( res.velocity.z, 0.0,     1e-2);
	}


	TEST(SimulationTest, ParticleMigrationOneParticle) {

		// this test checks whether particles are properly migrated between cells

		// Set universe properties
		UniverseProperties properties;
		properties.size = {2,1,1};
		properties.cellWidth = { .5,.5,.5 };
		properties.dt = 0.1;
		properties.useCase = UseCase::Test;

		// Create a universe with these properties
		Universe universe = Universe(properties);

		// configure the cell
		Cell& a = universe.cells[{0,0,0}];
		Cell& b = universe.cells[{1,0,0}];

		// initialize the field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.2, 0.0, 0.0 };
			field[pos].B = { 0.2, 0.0, 0.0 };
		});

		// add one particle
		Particle p;
		p.position.x = p.position.y = 0.4;
		p.position.z = 0;
		p.velocity.z = p.velocity.y = 0;
		p.velocity.x = 1;
		p.q = 1;
		p.mass = 1;

		// add test particle to first cell
		a.particles.push_back(p);

		// check particle position
		ASSERT_FALSE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());

		// number of steps
		unsigned niter = 4;

		// run the simulation
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check whether the particle was moved from one cell to another
		ASSERT_TRUE(a.particles.empty());
		ASSERT_FALSE(b.particles.empty());

		// run the simulation
		niter = 1;
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check whether the particle was moved out of the domain
		// TODO: the particle has to move back to cell a due to the periodic boundary conditions
		ASSERT_FALSE(a.particles.empty());
		ASSERT_TRUE(b.particles.empty());
	}


	TEST(SimulationTest, ParticleMigrationManyParticles) {

		// this test checks whether particles are properly migrated between cells

		// Set universe properties
		UniverseProperties properties;
		properties.size = {2,2,2};
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
		p0.velocity.z = 0.6; p0.velocity.y = p0.velocity.z = 0.6;
		p0.q = p0.mass = 1;
		p1.position.x = 0.8; p1.position.y = p1.position.z = 0.4;
		p1.velocity.z = 0.3; p1.velocity.y = p1.velocity.z = 0.3;
		p1.q = p1.mass = 1;
		p2.position.y = 0.8; p2.position.x = p2.position.z = 0.4;
		p2.velocity.z = 0.1; p2.velocity.y = p2.velocity.z = 0.1;
		p2.q = p2.mass = 1;
		p3.position.x = 0.75; p3.position.y = 0.75; p3.position.z = 0.4;
		p3.velocity.z = 0.8; p3.velocity.y = p3.velocity.z = 0.8;
		p3.q = p3.mass = 1;
		p4 = p0;
		p4.position.z = 0.6;
		p5 = p1;
		p5.position.z = 0.8;
		p6 = p3;
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

		// number of steps
		unsigned niter = 20;

		// run the simulation
		simulateSteps<detail::default_particle_to_field_projector, detail::default_field_solver, detail::boris_mover>(niter,universe);

		// check number of particles in the domain
		std::atomic<int> total_particles = ATOMIC_VAR_INIT(0);;
		decltype(field.size()) size = 2;
		allscale::api::user::pfor(zero,size,[&](auto& pos){
			std::atomic_fetch_add( &total_particles, (int) universe.cells[pos].particles.size() );
		});
		EXPECT_EQ(8, total_particles);
	}

} // end namespace ipic3d
