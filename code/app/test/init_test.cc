#include <gtest/gtest.h>

#include "ipic3d/app/init_properties.h"

namespace ipic3d {

	TEST(InitProperties, Basic) {

		using Vu = Vector3<unsigned>;
		using Vd = Vector3<double>;

	    // Set some universe properties
		InitProperties properties;

		EXPECT_EQ(1, properties.numSteps);
		EXPECT_EQ(std::vector<Vu>{}, properties.particlesPerCell);
		EXPECT_EQ(std::vector<Vd>{}, properties.driftVelocity);
		EXPECT_EQ((Vd{0.0,0.0,0.0}), properties.magneticField);
		EXPECT_EQ(1.0, properties.rhoInit);
		

		InitProperties init(100, {10,20,30}, {0.2,0.3,0.4}, {1.3,1.6,0.0001}, 1.23456);
		EXPECT_EQ(100, init.numSteps);
		EXPECT_EQ((std::vector<Vu>{10,20,30}), init.particlesPerCell);
		EXPECT_EQ((std::vector<Vd>{0.2,0.3,0.4}), init.driftVelocity);
		EXPECT_EQ((Vd{1.3,1.6,0.0001}), init.magneticField);
		EXPECT_EQ(1.23456, init.rhoInit);
    }

	TEST(InitProperties, Printable) {
		InitProperties properties;
		// ensures that InitProperties can be printed
		std::stringstream ss;
		ss << properties;
		EXPECT_FALSE(ss.str().empty());
	}

} // end namespace ipic3d
