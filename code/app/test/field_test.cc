#include <gtest/gtest.h>

#include "ipic3d/app/cell.h"
#include "ipic3d/app/universe.h"

namespace ipic3d {

	TEST(Field, curlB) {
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 1, 1, 1 };
		properties.cellWidth = { 1.0, 1.0, 1.0 };
		properties.useCase = UseCase::Dipole;

		utils::Coordinate<3> pos{1, 1, 1};

		// initialize field
		BcField bcfields(properties.size + coordinate_type(1));
		decltype(bcfields.size()) zero = 0;
		allscale::api::user::pfor(zero,bcfields.size(),[&](auto& pos){
			bcfields[pos].Bc = { 0.0, 0.0, 0.0 };
		});

		Vector3<double> curlB = { 0.0, 0.0, 0.0 };

		computeCurlB(properties, pos, bcfields, curlB); 
		EXPECT_NEAR( curlB.x, 0.0, 1e-06 );
		EXPECT_NEAR( curlB.y, 0.0, 1e-06 );
		EXPECT_NEAR( curlB.z, 0.0, 1e-06 );

		// change the values of electric field
	    for(int i = 0; i < 2; i++) {
		    for(int j = 0; j < 2; j++) {
			    for(int k = 0; k < 2; k++) {
					utils::Coordinate<3> cur({pos[0]-i,pos[1]-j,pos[2]-k});
				    bcfields[cur].Bc = {0.0, 2.0, i * 1.0 + j * 2.0 + k * 3.0};
			    }
		    }
	    }

		// re-evaluate
		computeCurlB(properties, pos, bcfields, curlB); 
		EXPECT_NEAR( curlB.x, -2.0, 1e-06 );
		EXPECT_NEAR( curlB.y, 1.0, 1e-06 );
		EXPECT_NEAR( curlB.z, 0.0, 1e-06 );

		// change the width of cells
		properties.cellWidth = { 1.0, 5.0, 10.0 };

		// re-evaluate
		computeCurlB(properties, pos, bcfields, curlB); 
		EXPECT_NEAR( curlB.x, -0.4, 1e-06 );
		EXPECT_NEAR( curlB.y, 1.0, 1e-06 );
		EXPECT_NEAR( curlB.z, 0.0, 1e-06 );
    }

	TEST(Field, curlE) {
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 1, 1, 1 };
		properties.cellWidth = { 1.0, 1.0, 1.0 };
		properties.useCase = UseCase::Dipole;

		utils::Coordinate<3> pos{0, 0, 0};

		// initialize field
		Field fields(properties.size + coordinate_type(1));
		decltype(fields.size()) zero = 0;
		allscale::api::user::pfor(zero,fields.size(),[&](auto& pos){
			fields[pos].E = { 0.0, 0.0, 0.0 };
		});

		Vector3<double> curlE = { 0.0, 0.0, 0.0 };

		computeCurlE(properties, pos, fields, curlE); 
		EXPECT_NEAR( curlE.x, 0.0, 1e-06 );
		EXPECT_NEAR( curlE.y, 0.0, 1e-06 );
		EXPECT_NEAR( curlE.z, 0.0, 1e-06 );

		// change the values of electric field
	    for(int i = 0; i < 2; i++) {
		    for(int j = 0; j < 2; j++) {
			    for(int k = 0; k < 2; k++) {
					utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});
				    fields[cur].E = {0.0, 2.0, i * 1.0 + j * 2.0 + k * 3.0};
			    }
		    }
	    }

		// re-evaluate
		computeCurlE(properties, pos, fields, curlE); 
		EXPECT_NEAR( curlE.x, 2.0, 1e-06 );
		EXPECT_NEAR( curlE.y, -1.0, 1e-06 );
		EXPECT_NEAR( curlE.z, 0.0, 1e-06 );

		// change the width of cells
		properties.cellWidth = { 1.0, 5.0, 10.0 };

		// re-evaluate
		computeCurlE(properties, pos, fields, curlE); 
		EXPECT_NEAR( curlE.x, 0.4, 1e-06 );
		EXPECT_NEAR( curlE.y, -1.0, 1e-06 );
		EXPECT_NEAR( curlE.z, 0.0, 1e-06 );
    }

	TEST(Field, interpC2N) {
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 1, 1, 1 };

	    // use center position
		utils::Coordinate<3> pos{1, 1, 1};

		// initialize field
		BcField bcfields(properties.size + coordinate_type(1));
		decltype(bcfields.size()) zero = 0;
		allscale::api::user::pfor(zero,bcfields.size(),[&](auto& pos){
			bcfields[pos].Bc = { 1.0, 1.0, 1.0 };
		});

		Field fields(properties.size);
		fields[pos].B = { 0.0, 0.0, 0.0};

		
		interpC2N(pos, bcfields, fields); 
		auto B = fields[pos].B;
		EXPECT_NEAR( B.x, 1.0, 1e-06 );
		EXPECT_NEAR( B.y, 1.0, 1e-06 );
		EXPECT_NEAR( B.z, 1.0, 1e-06 );

		// change the values of magnetic field on nodes
	    for(int i = 0; i < 2; i++) {
		    for(int j = 0; j < 2; j++) {
			    for(int k = 0; k < 2; k++) {
					utils::Coordinate<3> cur({pos[0]-i,pos[1]-j,pos[2]-k});
				    bcfields[cur].Bc = {0.0, (i + 1.0) * (j + 1.0) * (k + 1.0), i * 1.0 + j * 2.0 + k * 3.0};
			    }
		    }
	    }

		// re-evaluate
		interpC2N(pos, bcfields, fields); 
		B = fields[pos].B;
		EXPECT_NEAR( B.x, 0.0, 1e-06 );
		EXPECT_NEAR( B.y, 3.375, 1e-06 );
		EXPECT_NEAR( B.z, 3.0, 1e-06 );
    }

	TEST(Field, interpN2C) {
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 1, 1, 1 };

	    // use center position
		utils::Coordinate<3> pos{0, 0, 0};

		// initialize field
		Field fields(properties.size + coordinate_type(1));
		decltype(fields.size()) zero = 0;
		allscale::api::user::pfor(zero,fields.size(),[&](auto& pos){
			fields[pos].B = { 1.0, 1.0, 1.0 };
		});

		BcField bcfields(properties.size);
		bcfields[pos].Bc = { 0.0, 0.0, 0.0};

		interpN2C(pos, fields, bcfields); 
		auto Bc = bcfields[pos].Bc;
		EXPECT_NEAR( Bc.x, 1.0, 1e-06 );
		EXPECT_NEAR( Bc.y, 1.0, 1e-06 );
		EXPECT_NEAR( Bc.z, 1.0, 1e-06 );

		// change the values of magnetic field on nodes
	    for(int i = 0; i < 2; i++) {
		    for(int j = 0; j < 2; j++) {
			    for(int k = 0; k < 2; k++) {
					utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});
				    fields[cur].B = {0.0, 2.0, i * 1.0 + j * 2.0 + k * 3.0};
			    }
		    }
	    }

		// re-evaluate
		interpN2C(pos, fields, bcfields); 
		Bc = bcfields[pos].Bc;
		EXPECT_NEAR( Bc.x, 0.0, 1e-06 );
		EXPECT_NEAR( Bc.y, 2.0, 1e-06 );
		EXPECT_NEAR( Bc.z, 3.0, 1e-06 );
    }

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
		Field field2(properties.size + coordinate_type(1));
		solveFieldStatically(properties, pos, field2);
		E = field2[pos].E;
		EXPECT_NEAR( E.x, 0.0, 1e-06 );
		EXPECT_NEAR( E.y, 0.0, 1e-06 );
		EXPECT_NEAR( E.z, 0.0, 1e-06 );

		B = field2[pos].B;
		EXPECT_NEAR( B.x, -1545900104359.654, 1e-02 );
		EXPECT_NEAR( B.y, -809757197521.7235, 1e-02 );
		EXPECT_NEAR( B.z,  4449574903553.713, 1e-02 );


        // verify test case
		properties.useCase = UseCase::Test;
		properties.size = { 11,11,21 };
	    pos = {10, 10, 20};
		Field field3(properties.size + coordinate_type(1));
		allscale::api::user::pfor(zero,field3.size(),[&](auto& pos){
			field3[pos].E = { 0.0, 0.0, 0.0 };
			field3[pos].B = { 0.0, 0.0, 0.0 };
		});
		solveFieldStatically(properties, pos, field3);

		E = field3[pos].E;
		EXPECT_NEAR( E.x, 0.0, 1e-06 );
		EXPECT_NEAR( E.y, 0.0, 1e-06 );
		EXPECT_NEAR( E.z, 0.0, 1e-06 );

		B = field3[pos].B;
		EXPECT_NEAR( B.x, 0.0, 1e-06 );
		EXPECT_NEAR( B.y, 0.0, 1e-06 );
		EXPECT_NEAR( B.z, 0.0, 1e-06 );
    }

} // end namespace ipic3d

