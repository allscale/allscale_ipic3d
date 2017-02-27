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
			p.velocity.x = 1;
			p.velocity.y = p.velocity.z = 0;
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

			EXPECT_EQ(1.0, a.particles.front().velocity.x);
			EXPECT_EQ(0.0, a.particles.front().velocity.y);
			EXPECT_EQ(0.0, a.particles.front().velocity.z);

			EXPECT_EQ(1.0, a.particles.front().mass);
			EXPECT_EQ(1.0, a.particles.front().q);

			// run one simulation step
			UseCase test = UseCase::Test;
			simulateStep<detail::default_particle_to_field_projector,detail::default_field_solver,Mover>(test,universe);

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

			universe.properties.dt = 0.25;
			int niter = 4;

			simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,Mover>(test,niter,universe);

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

	TEST(SimulationTest, SingleParticleFirstOrder) {

		// this test checks whether particles are properly migrated between cells
		detail::testSingleParticle<detail::default_particle_mover>();

	}

	TEST(SimulationTest, SingleParticleSecondOrder) {

		// this test checks whether particles are properly migrated between cells
		detail::testSingleParticle<detail::boris_mover>();

	}


	TEST(SimulationTest, SingleParticleBorisMover) {

		// Set universe properties
		UniverseProperties properties;
		properties.size = {1,1,1};
		properties.cellWidth = { 100,100,100 };
		properties.dt = 0.1;

		// number of steps
		int niter = 9;

		// Create a universe with these properties
		Universe universe = Universe(properties);

		// configure the cell
		Cell& cell = universe.cells[{0,0,0}];

		// initialize the field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.2, 0, 0 };
			field[pos].B = { 0.2, 0, 0 };
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
		UseCase test = UseCase::Test;
		simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(test,niter,universe);

		// check where particle ended up
		ASSERT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected
		EXPECT_NEAR( res.position.x, 0.590, 0.001);
		EXPECT_NEAR( res.position.y, 0.589, 0.001);
		EXPECT_NEAR( res.position.z, 0.894, 0.001);

	}


	TEST(SimulationTest, SingleParticleBorisMoverLarmorRadiusPseudo) {

		// Set universe properties
		UniverseProperties properties;
		properties.size = { 1,2,1 };
		properties.cellWidth = { 10,5,10 };
		properties.dt = 0.1;

		// number of steps
		int niter = 100;

		// Create Universe with these properties
		Universe universe = Universe(properties);

		Cell& a = universe.cells[{0,0,0}];
		Cell& b = universe.cells[{0,1,0}];

		// initialize field
		Field& field = universe.field;
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.0, 0.0, 0.0 };
			field[pos].B = { 0.0, 0.0, 0.1 };
		});

		// add one particle
		Particle p;
		p.position.x = p.position.y = p.position.z = 0.0;

		p.velocity.x = p.velocity.z = 0.0;
		p.velocity.y = 1.0;

		p.q = -1.0;
		p.mass = 1.0;

		// compute Larmor radius
		double rL = p.mass * p.velocity.y / (fabs(p.q) * field[{0,0,0}].B.z);
		EXPECT_NEAR( rL, 10, 0.1 );
		p.position.x = rL;

		// push velocity back in time by 1/2 dt
		p.updateVelocityBorisStyle(field[{0,0,0}].E, field[{0,0,0}].B, -0.5*properties.dt);

		a.particles.push_back(p);

		// check particle position
		EXPECT_FALSE(a.particles.empty());
		EXPECT_TRUE(b.particles.empty());

		// run the simulation
		UseCase test = UseCase::Test;
		simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(test,niter,universe);

		// check where particle ended up
		EXPECT_TRUE(a.particles.empty());
		EXPECT_FALSE(b.particles.empty());

		Particle res = b.particles.front();

		// check that the position is close to what is expected
		// comparing against the matlab code after 10 iterations
		EXPECT_NEAR( res.position.x, 5.4030, 1e-04);
		EXPECT_NEAR( res.position.y, 8.4147, 1e-04);
		EXPECT_NEAR( res.position.z, 0.0,	 1e-05);

		// check that the velocity is close to what is expected
		EXPECT_NEAR( res.velocity.x, -0.8388, 1e-04);
		EXPECT_NEAR( res.velocity.y,  0.5445, 1e-04);
		EXPECT_NEAR( res.velocity.z,  0.0,    1e-04);
	}


	TEST(SimulationTest, SingleParticleBorisMoverLarmorRadius) {

		// Set universe properties
		UniverseProperties properties;
		properties.size = { 1,1,1 };
		properties.cellWidth = { 1e4,1e4,1e4 };
		properties.dt = 3e-11;

		// number of steps
		int niter = 10;

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
		UseCase test = UseCase::Test;
		simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(test,niter,universe);

		// check where particle ended up
		EXPECT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected
		// comparing against the matlab code after 10 iterations
		EXPECT_NEAR( res.position.x, 4.9127e-05, 1e-06);
		EXPECT_NEAR( res.position.y, 2.8631e-05, 1e-06);
		EXPECT_NEAR( res.position.z, 0.0,		 1e-06);

		// check that the velocity is close to what is expected
		EXPECT_NEAR( res.velocity.x, -4.8050e+04, 1e2);
		EXPECT_NEAR( res.velocity.y,  8.7699e+04, 1e2);
		EXPECT_NEAR( res.velocity.z,  0.0,        1e-2);
	}

} // end namespace ipic3d
