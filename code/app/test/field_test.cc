#include <gtest/gtest.h>

#include "ipic3d/app/cell.h"

namespace ipic3d {

	TEST(FieldTest, TrilinearInterpolation) {
	    // use center position
	    Vector3<double> pos{0.5, 0.5, 0.5};

	    // create 8 points, all with the same value
	    double fieldCorners[2][2][2];
	    for(int i = 0; i < 2; i++) {
		    for(int j = 0; j < 2; j++) {
			    for(int k = 0; k < 2; k++) {
				    fieldCorners[i][j][k] = 1.0;
			    }
		    }
	    }

	    // perform trilinear interpolation and check result
	    double fieldAtPos = trilinearInterpolation(fieldCorners, pos);
	    EXPECT_EQ(fieldAtPos, 1.0);

	    // change values of points
	    for(int i = 0; i < 2; i++) {
		    for(int j = 0; j < 2; j++) {
			    for(int k = 0; k < 2; k++) {
				    fieldCorners[i][j][k] = i * 1 + j * 2 + k * 3;
			    }
		    }
	    }

	    // re-evaluate
	    fieldAtPos = trilinearInterpolation(fieldCorners, pos);
	    EXPECT_EQ(fieldAtPos, 3.0);

	    // change target position and re-evaluate
	    pos = {0.0, 0.0, 0.0};
	    fieldAtPos = trilinearInterpolation(fieldCorners, pos);
	    EXPECT_EQ(fieldAtPos, 0.0);
	    pos = {1.0, 1.0, 1.0};
	    fieldAtPos = trilinearInterpolation(fieldCorners, pos);
	    EXPECT_EQ(fieldAtPos, 6.0);
    }

} // end namespace ipic3d

