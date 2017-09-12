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
		p.updatePosition(1);

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
		p.updateVelocity({ 0.0, 0.0, 0.0 }, 1);

		EXPECT_EQ(1.0, p.velocity.x);
		EXPECT_EQ(2.0, p.velocity.y);
		EXPECT_EQ(3.0, p.velocity.z);

		// apply some force
		p.updateVelocity({ 1.0, -1.0, 3.0 }, 1);

		EXPECT_EQ(1.0 + ( 1.0 * p.qom), p.velocity.x);
		EXPECT_EQ(2.0 + (-1.0 * p.qom), p.velocity.y);
		EXPECT_EQ(3.0 + ( 3.0 * p.qom), p.velocity.z);


	}

	TEST(Particle, Force) {

		// this test checks proper force calculation of a particle in an electric field

		// create a test-particle
		Particle p;
		p.position.x = 0.5;
		p.position.y = 0.5;
		p.position.z = 0.5;

		p.q = 0.25;

		p.velocity.x = 0.0;
		p.velocity.y = 0.0;
		p.velocity.z = 0.0;

		// initialize an electric field
		Vector3<double> E[2][2][2];
		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 2; j++) {
				for(int k = 0; k < 2; k++) {
					E[i][j][k] = 0;
				}
			}
		}

		Vector3<double> force = computeElectricForce(E, p);

		EXPECT_EQ(0.0, force.x);
		EXPECT_EQ(0.0, force.y);
		EXPECT_EQ(0.0, force.z);

		// modify electric field
		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 2; j++) {
				for(int k = 0; k < 2; k++) {
					E[i][j][k] = 1.0;
				}
			}
		}

		force = computeElectricForce(E, p);

		EXPECT_EQ(2.0, force.x);
		EXPECT_EQ(2.0, force.y);
		EXPECT_EQ(2.0, force.z);

	}


} // end namespace ipic3d

