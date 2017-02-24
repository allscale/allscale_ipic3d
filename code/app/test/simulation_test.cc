#include <gtest/gtest.h>

#include "ipic3d/app/simulator.h"
#include "ipic3d/app/universe.h"

#include "allscale/api/user/operator/pfor.h"

namespace ipic3d {


	namespace detail {


		template<typename Mover>
		void testSingleParticle() {

			// this test checks whether particles are properly migrated between cells

			// Create Universe with two cells and corresponding Field
			Universe universe = Universe({ 2,1,1 });

			Field& field = universe.field;
			decltype(field.size()) zero = 0;
			allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
				field[pos].E = { 0, 0, 0 };
				field[pos].B = { 0, 0, 0 };
			});

			Cell& a = universe.cells[{0,0,0}];
			Cell& b = universe.cells[{1,0,0}];

			// fix cell positions
			a.center.x = 0.5;
			b.center.x = 1.5;
			a.center.y = a.center.z = b.center.y = b.center.z = 0.5;

			// fix cell widths
			a.spacing.x = b.spacing.x = 1;
			a.spacing.y = b.spacing.y = 1;
			a.spacing.z = b.spacing.z = 1;

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
			simulateStep<detail::default_particle_to_field_projector,detail::default_field_solver,Mover>(1,universe);

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

			simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,Mover>(4,0.25,universe);

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

		// Create Universe with one Cell and corresponding Field
		Universe universe = Universe({ 1,1,1 });

		// configure the cell
		Cell& cell = universe.cells[{0,0,0}];
		cell.center.x = cell.center.y = cell.center.z = 50;
		cell.spacing.x = cell.spacing.y = cell.spacing.z = 100;

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
		simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(9,0.1,universe);

		// check where particle ended up
		ASSERT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected
		EXPECT_NEAR( res.position.x, 0.590, 0.001);
		EXPECT_NEAR( res.position.y, 0.589, 0.001);
		EXPECT_NEAR( res.position.z, 0.894, 0.001);

	}


	TEST(SimulationTest, SingleParticleBorisMoverLarmorRadiusPseudo) {
		int niter = 100;
		double dt = 0.1;

		// Create Universe with one Cell and corresponding Field
		Universe universe = Universe({1,2,1});

		Cell& a = universe.cells[{0,0,0}];
		Cell& b = universe.cells[{0,1,0}];

		// configure cells
		a.center.y = 2.5;
		b.center.y = 7.5;
		a.center.x = a.center.z = b.center.x = b.center.z = 5.0;

		// fix cell widths
		a.spacing.y = b.spacing.y =  5.0;
		a.spacing.x = b.spacing.x = 10.0;
		a.spacing.z = b.spacing.z = 10.0;

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
		p.updateVelocityBorisStyle(field[{0,0,0}].E, field[{0,0,0}].B, -0.5*dt);

		a.particles.push_back(p);

		// check particle position
		EXPECT_FALSE(a.particles.empty());
		EXPECT_TRUE(b.particles.empty());

		// run the simulation
		simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(niter,dt,universe);

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
		int niter = 10;
		double dt = 3e-11;

		// Create Universe with one Cell and corresponding Field
		Universe universe = Universe({ 1,1,1 });
		Cell& cell = universe.cells[{0, 0, 0}];

		// configure the cell
		cell.center.x = cell.center.y = cell.center.z = 0;
		cell.spacing.x = cell.spacing.y = cell.spacing.z = 1e4;


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
		p.updateVelocityBorisStyle(field[{0,0,0}].E, field[{0,0,0}].B, -0.5*dt);

		cell.particles.push_back(p);

		// check particle position
		EXPECT_FALSE(cell.particles.empty());

		// run the simulation
		simulateSteps<detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(niter,dt,universe);

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
