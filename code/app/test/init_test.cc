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
		EXPECT_EQ((Vd{0,0,0}), properties.magneticFieldAmplitude);
    }

	TEST(InitProperties, Printable) {
		InitProperties properties;
		// ensures that InitProperties can be printed
		std::stringstream ss;
		ss << properties;
		EXPECT_FALSE(ss.str().empty());
	}

} // end namespace ipic3d
