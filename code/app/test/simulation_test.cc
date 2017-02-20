#include <gtest/gtest.h>

#include "ipic3d/app/simulator.h"

#include "allscale/api/user/operator/pfor.h"

namespace ipic3d {


	namespace detail {


		template<typename Mover>
		void testSingleParticle() {

			// this test checks whether particles are properly migrated between cells

			// create a grid of two cells
			Cells cells = Cells({2,1,1});
			Field fields = Field({3,2,2});

			decltype(fields.size()) zero = 0;
			allscale::api::user::pfor(zero,fields.size(),[&](auto& pos){
				fields[pos].E = { 0, 0, 0 };
				fields[pos].B = { 0, 0, 0 };
			});

			Cell& a = cells[{0,0,0}];
			Cell& b = cells[{1,0,0}];

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
			EXPECT_FALSE(a.particles.empty());
			EXPECT_TRUE(b.particles.empty());

			// run one simulation step
			simulateStep<Cell,FieldNode,detail::default_particle_to_field_projector,detail::default_field_solver,Mover>(1,cells,fields);

			EXPECT_TRUE(a.particles.empty());
			EXPECT_FALSE(b.particles.empty());

			// -- check the particle position
			EXPECT_EQ(1.5,b.particles.front().position.x);
			EXPECT_EQ(0.5,b.particles.front().position.y);
			EXPECT_EQ(0.5,b.particles.front().position.z);

			EXPECT_EQ(1.0,b.particles.front().velocity.x);
			EXPECT_EQ(0.0,b.particles.front().velocity.y);
			EXPECT_EQ(0.0,b.particles.front().velocity.z);

			// change velocity and send back
			b.particles.front().velocity.x = -1;

			simulateSteps<Cell,FieldNode,detail::default_particle_to_field_projector,detail::default_field_solver,Mover>(4,0.25,cells,fields);

			EXPECT_FALSE(a.particles.empty());
			EXPECT_TRUE(b.particles.empty());

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

		// create one cell
		Cells cells({1,1,1});
		Cell& cell = cells[{0,0,0}];

		// configure the cell
		cell.center.x = cell.center.y = cell.center.z = 50;
		cell.spacing.x = cell.spacing.y = cell.spacing.z = 100;


		// create a surrounding force field
		Field fields({2,2,2});

		decltype(fields.size()) zero = 0;
		allscale::api::user::pfor(zero,fields.size(),[&](auto& pos){
			fields[pos].E = { 0.2, 0, 0 };
			fields[pos].B = { 0.2, 0, 0 };
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
		simulateSteps<Cell,FieldNode,detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(9,0.1,cells,fields);

		// check where particle ended up
		ASSERT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected
		EXPECT_NEAR( res.position.x, 0.590, 0.001);
		EXPECT_NEAR( res.position.y, 0.589, 0.001);
		EXPECT_NEAR( res.position.z, 0.894, 0.001);

	}


	TEST(SimulationTest, SingleParticleBorisMoverLarmorRadiusPseudo) {
		int niter = 10;
		double dt = 0.1;

		// create one cell
		Cells cells({1,1,1});
		Cell& cell = cells[{0,0,0}];

		// configure the cell
		cell.center.x = cell.center.y = cell.center.z = 5;
		cell.spacing.x = cell.spacing.y = cell.spacing.z = 10;


		// create a surrounding force field
		Field fields({2,2,2});

		decltype(fields.size()) zero = 0;
		allscale::api::user::pfor(zero,fields.size(),[&](auto& pos){
			fields[pos].E = { 0.0, 0.0, 0.0 };
			fields[pos].B = { 0.0, 0.0, 0.1 };
		});

		// add one particle
		Particle p;
		p.position.x = p.position.y = p.position.z = 0.0;

		p.velocity.x = p.velocity.z = 0.0;
		p.velocity.y = 1;

		p.q = -1;
		p.mass = 1;

		// compute Larmor radius
		double rL = p.mass * p.velocity.y / (fabs(p.q) * fields[{0,0,0}].B.z);
		EXPECT_NEAR( rL, 10, 0.1 );
		p.position.x = rL;

		// push velocity back in time by 1/2 dt
		p.updateVelocityBorisStyle(fields[{0,0,0}].E, fields[{0,0,0}].B, -0.5*dt);

		cell.particles.push_back(p);

		// run the simulation
		simulateSteps<Cell,FieldNode,detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(niter,dt,cells,fields);

		// check where particle ended up
		ASSERT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected 
		// comparing against the matlab code after 10 iterations
		EXPECT_NEAR( res.position.x, 9.9500, 0.0001);
		EXPECT_NEAR( res.position.y, 0.9983, 0.0001);
		EXPECT_NEAR( res.position.z, 0.0,	0.00001);

		// check that the velocity is close to what is expected
		EXPECT_NEAR( res.velocity.x, -0.0948, 0.0001);
		EXPECT_NEAR( res.velocity.y, 0.9954,  0.0001);
		EXPECT_NEAR( res.velocity.z, 0.0,     0.0001);
	}


	TEST(SimulationTest, SingleParticleBorisMoverLarmorRadius) {
		int niter = 10;
		double dt = 0.1;

		// create one cell
		Cells cells({1,1,1});
		Cell& cell = cells[{0,0,0}];

		// configure the cell
		cell.center.x = cell.center.y = cell.center.z = 5;
		cell.spacing.x = cell.spacing.y = cell.spacing.z = 10;


		// create a surrounding force field
		Field fields({2,2,2});

		decltype(fields.size()) zero = 0;
		allscale::api::user::pfor(zero,fields.size(),[&](auto& pos){
			fields[pos].E = { 0.0, 0.0, 0.0  };
			fields[pos].B = { 0.0, 0.0, 0.01 };
		});

		// add one particle
		Particle p;
		p.position.x = p.position.y = p.position.z = 0.0;

		p.velocity.x = p.velocity.z = 0.0;
		p.velocity.y = 1e5;

		p.q = -1.602e-19;
		p.mass = 9.109e-31;

		// compute Larmor radius
		double rL = p.mass * p.velocity.y / (fabs(p.q) * fields[{0,0,0}].B.z);
		EXPECT_NEAR( rL, 5.686e-05, 1e-06 );
		p.position.x = rL;
	
		return;
		// push velocity back in time by 1/2 dt
		p.updateVelocityBorisStyle(fields[{0,0,0}].E, fields[{0,0,0}].B, -0.5*dt);

		cell.particles.push_back(p);

		// run the simulation
		simulateSteps<Cell,FieldNode,detail::default_particle_to_field_projector,detail::default_field_solver,detail::boris_mover>(niter,dt,cells,fields);

		// check where particle ended up
		ASSERT_FALSE(cell.particles.empty());
		Particle res = cell.particles.front();

		// check that the position is close to what is expected 
		// comparing against the matlab code after 10 iterations
		EXPECT_NEAR( res.position.x, 9.9500, 0.0001);
		EXPECT_NEAR( res.position.y, 0.9983, 0.0001);
		EXPECT_NEAR( res.position.z, 0.0,	0.00001);

		// check that the velocity is close to what is expected
		EXPECT_NEAR( res.velocity.x, -0.0948, 0.0001);
		EXPECT_NEAR( res.velocity.y, 0.9954,  0.0001);
		EXPECT_NEAR( res.velocity.z, 0.0,     0.0001);
	}

} // end namespace ipic3d
