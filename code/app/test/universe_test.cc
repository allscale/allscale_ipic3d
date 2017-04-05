#include <gtest/gtest.h>

#include "ipic3d/app/universe.h"

namespace ipic3d {

	TEST(UniverseProperties, Basic) {

	    using V = Vector3<double>;

	    // Set some universe properties
		UniverseProperties properties;

		EXPECT_EQ(UseCase::ParticleWave, properties.useCase);
	    EXPECT_EQ((coordinate_type{1, 1, 1}), properties.size);
	    EXPECT_EQ((V{1, 1, 1}), properties.cellWidth);
	    EXPECT_EQ(1, properties.dt);
	    EXPECT_EQ(0.0, properties.objectRadius);
	    EXPECT_EQ((V{0, 0, 0}), properties.objectCenter);
	    EXPECT_EQ((V{0, 0, 0}), properties.magneticField);

		UniverseProperties properties2(UseCase::Test, { 2,3,4 }, { 0.2,0.3,0.4 }, 5, 42, { 0.7,0.8,0.9 }, {-0.1,-0.2,-0.3});

		EXPECT_EQ(UseCase::Test, properties2.useCase);
		EXPECT_EQ((coordinate_type{ 2,3,4 }), properties2.size);
		EXPECT_EQ((V{ 0.2,0.3,0.4 }), properties2.cellWidth);
		EXPECT_EQ(5, properties2.dt);
		EXPECT_EQ(42, properties2.objectRadius);
		EXPECT_EQ((V{ 0.7,0.8,0.9 }), properties2.objectCenter);
		EXPECT_EQ((V{ -0.1,-0.2,-0.3 }), properties2.magneticField);
    }

	TEST(UniverseProperties, Printable) {
		UniverseProperties properties;
		// ensures that UniverseProperties can be printed
		std::stringstream ss;
		ss << properties;
		EXPECT_FALSE(ss.str().empty());
	}

	TEST(UseCase, Printable) {
		UseCase useCase = UseCase::ParticleWave;
		std::stringstream ss;
		ss << useCase;
		EXPECT_EQ("ParticleWave", ss.str());
	}

	TEST(UniverseProperties, TypeProperties) {
		UniverseProperties properties;
		Universe universe(properties);
		EXPECT_TRUE(std::is_const<decltype(universe.properties)>::value);
	}

	TEST(Universe, Size) {

		// Set some universe properties
		UniverseProperties properties;
		properties.size = { 1,2,3 };

		// Create a universe with these properties
		Universe universe = Universe(properties);

		EXPECT_EQ(1, universe.properties.size.x);
		EXPECT_EQ(2, universe.properties.size.y);
		EXPECT_EQ(3, universe.properties.size.z);

		// Check the properties of the universe
		EXPECT_EQ(1, universe.properties.cellWidth.x);
		EXPECT_EQ(1, universe.properties.cellWidth.y);
		EXPECT_EQ(1, universe.properties.cellWidth.z);

		EXPECT_EQ(1, universe.properties.dt);

		// Get cells and field
		Cells& cells = universe.cells;
	    Field& field = universe.field;

		// Check its dimensions for cells and field
	    ASSERT_EQ(coordinate_type({1, 2, 3}), cells.size());
		// field grid should be 1 larger than cell grid in every dimension
	    ASSERT_EQ(coordinate_type({4, 5, 6}), field.size());

		// check a few cell positions
	    using Point = Vector3<double>;
	    EXPECT_EQ(Point({0.5, 0.5, 0.5}), getCenterOfCell({0, 0, 0}, properties));
	    EXPECT_EQ(Point({0.5, 0.5, 1.5}), getCenterOfCell({0, 0, 1}, properties));
	    EXPECT_EQ(Point({0.5, 1.5, 2.5}), getCenterOfCell({0, 1, 2}, properties));
    }

	TEST(Universe, Move) {

		EXPECT_TRUE(std::is_move_constructible<Universe>::value);

	}


} // end namespace ipic3d
