#include <gtest/gtest.h>

#include "ipic3d/app/universe.h"

namespace ipic3d {

	TEST(UniverseProperties, Basic) {

	    using V = Vector3<double>;

	    // Set some universe properties
		UniverseProperties properties;

		EXPECT_EQ(UseCase::Dipole, properties.useCase);
	    EXPECT_EQ((coordinate_type{1, 1, 1}), properties.size);
	    EXPECT_EQ((V{1.0, 1.0, 1.0}), properties.cellWidth);
	    EXPECT_EQ(1.0, properties.dt);
	    EXPECT_EQ(1.0, properties.speedOfLight);
	    EXPECT_EQ(0.0, properties.planetRadius);
	    EXPECT_EQ((V{0.0, 0.0, 0.0}), properties.objectCenter);
	    EXPECT_EQ((V{0.0, 0.0, 0.0}), properties.origin);
		EXPECT_EQ((V{0.0,0.0,0.0}), properties.externalMagneticField);
	    EXPECT_EQ(100, properties.FieldOutputCycle);
	    EXPECT_EQ(100, properties.ParticleOutputCycle);


		UniverseProperties properties2(UseCase::Dipole, { 3,3,3 }, { 0.2,0.3,0.4 }, 5.0, 13.0, 42.0, { 0.7,0.8,0.9 }, { -2.2, -2.5, -3.6 }, { 0.0,2.0,3.6 }, 100, 50);

		EXPECT_EQ(UseCase::Dipole, properties2.useCase);
		EXPECT_EQ((coordinate_type{ 3,3,3 }), properties2.size);
		EXPECT_EQ((V{ 0.2,0.3,0.4 }), properties2.cellWidth);
		EXPECT_EQ(5.0, properties2.dt);
		EXPECT_EQ(13.0, properties2.speedOfLight);
		EXPECT_EQ(42.0, properties2.planetRadius);
		EXPECT_EQ((V{ 0.7,0.8,0.9 }), properties2.objectCenter);
		EXPECT_EQ((V{ -2.2,-2.5,-3.6 }), properties2.origin);
		EXPECT_EQ((V{0.0,2.0,3.6}), properties2.externalMagneticField);
	    EXPECT_EQ(100, properties2.FieldOutputCycle);
	    EXPECT_EQ(50, properties2.ParticleOutputCycle);
    }

	TEST(UniverseProperties, Printable) {
		UniverseProperties properties;
		// ensures that UniverseProperties can be printed
		std::stringstream ss;
		ss << properties;
		EXPECT_FALSE(ss.str().empty());
	}

	TEST(UseCase, Printable) {
		UseCase useCase2 = UseCase::Dipole;
		std::stringstream ss2;
		ss2 << useCase2;
		EXPECT_EQ("Dipole", ss2.str());
	}

	TEST(UniverseProperties, TypeProperties) {
		UniverseProperties properties;
		Universe universe(properties);
		EXPECT_TRUE(std::is_const<decltype(universe.properties)>::value);
	}

	TEST(Universe, Size) {

		// Set some universe properties
		UniverseProperties properties;
		properties.size = { 4,4,4 };

		// Create a universe with these properties
		Universe universe = Universe(properties);

		EXPECT_EQ(4, universe.properties.size.x);
		EXPECT_EQ(4, universe.properties.size.y);
		EXPECT_EQ(4, universe.properties.size.z);

		// Check the properties of the universe
		EXPECT_EQ(1, universe.properties.cellWidth.x);
		EXPECT_EQ(1, universe.properties.cellWidth.y);
		EXPECT_EQ(1, universe.properties.cellWidth.z);

		EXPECT_EQ(1, universe.properties.dt);
		EXPECT_EQ(1, universe.properties.speedOfLight);

		// Get cells and field
		Cells& cells = universe.cells;
	    Field& field = universe.field;

		// Check its dimensions for cells and field
	    ASSERT_EQ(coordinate_type({4, 4, 4}), cells.size());
		// field grid should be 1 larger than cell grid in every dimension
	    ASSERT_EQ(coordinate_type({7, 7, 7}), field.size());

		// check a few cell positions
	    using Point = Vector3<double>;
	    EXPECT_EQ(Point({0.5, 0.5, 0.5}), getCenterOfCell({0, 0, 0}, properties));
	    EXPECT_EQ(Point({0.5, 0.5, 1.5}), getCenterOfCell({0, 0, 1}, properties));
	    EXPECT_EQ(Point({0.5, 1.5, 2.5}), getCenterOfCell({0, 1, 2}, properties));
    }

	TEST(Universe, Move) {

		EXPECT_TRUE(std::is_move_constructible<Universe>::value);

	}

	TEST(Universe, createUniverseFromParams) {

		// this test checks the createUniverseFromParams() function

		std::string path = std::string(PATH_TO_INPUTS) + "/test.inp";
		auto params = Parameters(path);

		// initialize initial properties
		Universe universe = createUniverseFromParams(params, "test");

		// verify the number of particles per cell
		int particlesPerCell = params.npcelx[0] * params.npcely[0] * params.npcelz[0];

		utils::Coordinate<3> zero = 0;
		allscale::api::user::algorithm::pfor(zero, universe.properties.size, [&](const utils::Coordinate<3>& pos) {
			EXPECT_EQ(particlesPerCell, (int) universe.cells[pos].particles.size());
		});
	}


} // end namespace ipic3d
