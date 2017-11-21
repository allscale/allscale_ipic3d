#include <gtest/gtest.h>

#include "ipic3d/app/particle.h"

namespace ipic3d {

	TEST(Particle, PositionUpdate) {

		// this test checks whether a particle is properly moved in space

		// create a test-particle
		Particle p;
		p.position.x = 0.0;
		p.position.y = 0.0;
		p.position.z = 0.0;

		p.velocity.x = 1.0;
		p.velocity.y = 2.0;
		p.velocity.z = 3.0;


		// check the particles position
		EXPECT_EQ(0.0, p.position.x);
		EXPECT_EQ(0.0, p.position.y);
		EXPECT_EQ(0.0, p.position.z);

		// update the particle position
		p.updatePosition(1.0);

		EXPECT_EQ(1.0, p.position.x);
		EXPECT_EQ(2.0, p.position.y);
		EXPECT_EQ(3.0, p.position.z);


		// update the particle position
		p.updatePosition(0.5);

		EXPECT_EQ(1.5, p.position.x);
		EXPECT_EQ(3.0, p.position.y);
		EXPECT_EQ(4.5, p.position.z);

	}


	TEST(Particle, VelocityUpdate) {

		// this test checks whether a particle is properly accelerated

		// create a test-particle
		Particle p;
		p.position.x = 0.0;
		p.position.y = 0.0;
		p.position.z = 0.0;

		p.qom = 1.0 / 0.25;

		p.velocity.x = 1.0;
		p.velocity.y = 2.0;
		p.velocity.z = 3.0;


		// check the particles velocity
		EXPECT_EQ(1.0, p.velocity.x);
		EXPECT_EQ(2.0, p.velocity.y);
		EXPECT_EQ(3.0, p.velocity.z);

		// no force, no change in velocity
		p.updateVelocity({ 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, 1.0);

		EXPECT_EQ(1.0, p.velocity.x);
		EXPECT_EQ(2.0, p.velocity.y);
		EXPECT_EQ(3.0, p.velocity.z);

		// apply some force
		p.updateVelocity({ 1.0, -1.0, 3.0 }, { 0.0, 1.0, 2.0 }, 1.0);

		EXPECT_NEAR(p.velocity.x, -2.42857, 1e-04);
		EXPECT_NEAR(p.velocity.y, 3.71429, 1e-04);
		EXPECT_NEAR(p.velocity.z, 12.14286, 1e-04);

	}

} // end namespace ipic3d

