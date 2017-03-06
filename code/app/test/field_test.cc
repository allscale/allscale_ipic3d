#include <gtest/gtest.h>

#include "ipic3d/app/cell.h"
#include "ipic3d/app/universe.h"

namespace ipic3d {

	TEST(Field, TrilinearInterpolation) {
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
        double vol = 1.0;
	    double fieldAtPos = trilinearInterpolationF2P(fieldCorners, pos, vol);
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
        vol = 0.75;
	    fieldAtPos = trilinearInterpolationF2P(fieldCorners, pos, vol);
	    EXPECT_EQ(fieldAtPos, 4.0);

	    // change target position and re-evaluate
	    pos = {0.0, 0.0, 0.0};
        vol = 3.0;
	    fieldAtPos = trilinearInterpolationF2P(fieldCorners, pos, vol);
	    EXPECT_EQ(fieldAtPos, 0.0);
	    pos = {1.0, 1.0, 1.0};
        vol = 1.0;
	    fieldAtPos = trilinearInterpolationF2P(fieldCorners, pos, vol);
	    EXPECT_EQ(fieldAtPos, 6.0);
    }

	TEST(Field, FieldSolverStatic) {
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 1,1,1 };
		properties.cellWidth = { 1,1,1 };
		properties.useCase = UseCase::Dipole;

		utils::Coordinate<3> pos{0, 0, 0};

		// initialize field
		Field field(properties.size + coordinate_type(1));
		decltype(field.size()) zero = 0;
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.0, 0.0, 0.0 };
			field[pos].B = { 0.0, 0.0, 0.0 };
		});

	    // apply static field solver and check results
		solveFieldStatically(properties, pos, field);
		auto E = field[pos].E;
		EXPECT_NEAR( E.x, 0.0, 1e-06 );
		EXPECT_NEAR( E.y, 0.0, 1e-06 );
		EXPECT_NEAR( E.z, 0.0, 1e-06 );

		auto B = field[pos].B;
		EXPECT_NEAR( B.x, -1.226388334605022e+16, 1e+02 );
		EXPECT_NEAR( B.y, -1.226388334605022e+16, 1e+02 );
		EXPECT_NEAR( B.z, 0.0, 1e-06 );


	    // change target position and re-evaluate
		properties.size = { 11,6,2 };
	    pos = {10, 5, 1};
		solveFieldStatically(properties, pos, field);
		E = field[pos].E;
		EXPECT_NEAR( E.x, 0.0, 1e-06 );
		EXPECT_NEAR( E.y, 0.0, 1e-06 );
		EXPECT_NEAR( E.z, 0.0, 1e-06 );

		B = field[pos].B;
		EXPECT_NEAR( B.x, -1545900104359.654, 1e-02 );
		EXPECT_NEAR( B.y, -809757197521.7235, 1e-02 );
		EXPECT_NEAR( B.z,  4449574903553.713, 1e-02 );

		properties.useCase = UseCase::Test;

        // verify test case
		allscale::api::user::pfor(zero,field.size(),[&](auto& pos){
			field[pos].E = { 0.0, 0.0, 0.0 };
			field[pos].B = { 0.0, 0.0, 0.0 };
		});

		properties.size = { 11,11,21 };
	    pos = {10, 10, 20};
		solveFieldStatically(properties, pos, field);

		E = field[pos].E;
		EXPECT_NEAR( E.x, 0.0, 1e-06 );
		EXPECT_NEAR( E.y, 0.0, 1e-06 );
		EXPECT_NEAR( E.z, 0.0, 1e-06 );

		B = field[pos].B;
		EXPECT_NEAR( B.x, 0.0, 1e-06 );
		EXPECT_NEAR( B.y, 0.0, 1e-06 );
		EXPECT_NEAR( B.z, 0.0, 1e-06 );
    }

} // end namespace ipic3d

