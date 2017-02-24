#include <gtest/gtest.h>

#include "ipic3d/app/universe.h"

namespace ipic3d {

	TEST(UniverseTest, Size) {

		// Set some universe properties
		UniverseProperties properties({1,2,3});

		// Create a universe with these properties
		Universe universe = Universe(properties);

		EXPECT_EQ(1, universe.properties.size.x);
		EXPECT_EQ(2, universe.properties.size.y);
		EXPECT_EQ(3, universe.properties.size.z);

		// Check the properties of the universe
		EXPECT_EQ(1, universe.properties.cellWidth.x);
		EXPECT_EQ(1, universe.properties.cellWidth.y);
		EXPECT_EQ(1, universe.properties.cellWidth.z);

		// Get cells and field
		Cells& cells = universe.cells;
	    Field& field = universe.field;

		// Check its dimensions for cells and field
	    ASSERT_EQ(coordinate_type({1, 2, 3}), cells.size());
		// field grid should be 1 larger than cell grid in every dimension
	    ASSERT_EQ(coordinate_type({2, 3, 4}), field.size());

		// check a few cell positions
	    using Point = Vector3<double>;
	    EXPECT_EQ(Point({0.5, 0.5, 0.5}), getCenterOfCell({0, 0, 0}, properties));
	    EXPECT_EQ(Point({1.5, 0.5, 1.5}), getCenterOfCell({1, 0, 1}, properties));
	    EXPECT_EQ(Point({1.5, 2.5, 3.5}), getCenterOfCell({1, 2, 3}, properties));
    }


} // end namespace ipic3d
