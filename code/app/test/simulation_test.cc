#include <gtest/gtest.h>

#include "ipic3d/simulator.h"

#include "allscale/api/user/operator/pfor.h"

namespace ipic3d {

	TEST(SimulationTest, SingleParticle) {

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
		a.x = 0.5;
		b.x = 1.5;
		a.y = a.z = b.y = b.z = 0.5;

		// fix cell widths
		a.dx = b.dx = 1;
		a.dy = b.dy = 1;
		a.dz = b.dz = 1;

		// insert one particle
		Particle p;
		p.x = p.y = p.z = 0.5;
		p.dx = 1;
		p.dy = p.dz = 0;
		p.mass = 1;
		p.q = 1;

		// add test particle to first cell
		a.particles.push_back(p);

		// check particle position
		EXPECT_FALSE(a.particles.empty());
		EXPECT_TRUE(b.particles.empty());

		// run one simulation step
		simulateStep(1,cells,fields);

		EXPECT_TRUE(a.particles.empty());
		EXPECT_FALSE(b.particles.empty());

		// -- check the particle position
		EXPECT_EQ(1.5,b.particles.front().x);
		EXPECT_EQ(0.5,b.particles.front().y);
		EXPECT_EQ(0.5,b.particles.front().z);

		EXPECT_EQ(1.0,b.particles.front().dx);
		EXPECT_EQ(0.0,b.particles.front().dy);
		EXPECT_EQ(0.0,b.particles.front().dz);

		// change velocity and send back
		b.particles.front().dx = -1;

		simulateSteps(4,0.25,cells,fields);

		EXPECT_FALSE(a.particles.empty());
		EXPECT_TRUE(b.particles.empty());

		EXPECT_EQ(0.5,a.particles.front().x);
		EXPECT_EQ(0.5,a.particles.front().y);
		EXPECT_EQ(0.5,a.particles.front().z);

		EXPECT_EQ(-1.0,a.particles.front().dx);
		EXPECT_EQ( 0.0,a.particles.front().dy);
		EXPECT_EQ( 0.0,a.particles.front().dz);

	}

} // end namespace ipic3d
