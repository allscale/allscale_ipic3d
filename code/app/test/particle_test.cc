#include <gtest/gtest.h>

#include "ipic3d/particle.h"

namespace ipic3d {

	TEST(SimulationTest, PositionUpdate) {

		// this test checks whether a particle is properly moved in space

		// create a test-particle
		Particle p;
		p.x = 0.0;
		p.y = 0.0;
		p.z = 0.0;

		p.dx = 1.0;
		p.dy = 2.0;
		p.dz = 3.0;


		// check the particles position
		EXPECT_EQ(0.0, p.x);
		EXPECT_EQ(0.0, p.y);
		EXPECT_EQ(0.0, p.z);

		// update the particle position
		p.updatePosition(1);

		EXPECT_EQ(1.0, p.x);
		EXPECT_EQ(2.0, p.y);
		EXPECT_EQ(3.0, p.z);


		// update the particle position
		p.updatePosition(0.5);

		EXPECT_EQ(1.5, p.x);
		EXPECT_EQ(3.0, p.y);
		EXPECT_EQ(4.5, p.z);

	}


	TEST(SimulationTest, VelocityUpdate) {

		// this test checks whether a particle is properly accelerated

		// create a test-particle
		Particle p;
		p.x = 0.0;
		p.y = 0.0;
		p.z = 0.0;

		p.mass = 0.25;

		p.dx = 1.0;
		p.dy = 2.0;
		p.dz = 3.0;


		// check the particles velocity
		EXPECT_EQ(1.0, p.dx);
		EXPECT_EQ(2.0, p.dy);
		EXPECT_EQ(3.0, p.dz);

		// no force, no change in velocity
		p.updateVelocity({ 0.0, 0.0, 0.0 }, 1);

		EXPECT_EQ(1.0, p.dx);
		EXPECT_EQ(2.0, p.dy);
		EXPECT_EQ(3.0, p.dz);

		// apply some force
		p.updateVelocity({ 1.0, -1.0, 3.0 }, 1);

		EXPECT_EQ(1.0 + ( 1.0 * 1 / p.mass), p.dx);
		EXPECT_EQ(2.0 + (-1.0 * 1 / p.mass), p.dy);
		EXPECT_EQ(3.0 + ( 3.0 * 1 / p.mass), p.dz);


	}



	TEST(SimulationTest, ParticleMigration) {
		// this test checks whether particles are properly migrated between cells
	}

} // end namespace ipic3d
