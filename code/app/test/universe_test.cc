#include <gtest/gtest.h>

#include "ipic3d/app/universe.h"

namespace ipic3d {

	TEST(UniverseTest, Size) {
		// Create a small universe
		Universe universe = Universe({1,2,3});

		// Get cells and field
		Cells& cells = universe.cells;
	    Field& field = universe.field;

		// Check its dimensions for cells and field
	    EXPECT_EQ(Universe::coordinate_type({1, 2, 3}), cells.size());
		// field grid should be 1 larger than cell grid in every dimension
	    EXPECT_EQ(Universe::coordinate_type({2, 3, 4}), field.size());
    }


} // end namespace ipic3d
