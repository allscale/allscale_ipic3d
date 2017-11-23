#include <gtest/gtest.h>

#include "ipic3d/app/cell.h"
#include "ipic3d/app/universe.h"

namespace ipic3d {

	TEST(Field, initField) {

		// this test verifies the init function

		std::string path = std::string(PATH_TO_INPUTS) + "/test.inp";
		auto params = Parameters(path);

		// initialize initial properties
		InitProperties initProperties = InitProperties(params);

		// initialize universe properties
		UniverseProperties universeProperties = UniverseProperties(params);


		// verify dipole init
		// to cover the else statement
		universeProperties.useCase = UseCase::Dipole;
		universeProperties.planetRadius = 1000; 
		Field fields = initFields(initProperties, universeProperties);

		auto start = utils::Coordinate<3>(1);
		auto end = universeProperties.size + coordinate_type(2);
		allscale::api::user::algorithm::pfor(start, end, [&](const utils::Coordinate<3>& cur) {
			EXPECT_NEAR(fields[cur].Bext.x, 0.0, 1e-15);
			EXPECT_NEAR(fields[cur].Bext.y, 0.0, 1e-15);
			EXPECT_NEAR(fields[cur].Bext.z, 0.0, 1e-15);
		});
	} 


	TEST(Field, curlB) {
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 2, 2, 2 };
		properties.cellWidth = { 1.0, 1.0, 1.0 };
		properties.useCase = UseCase::Dipole;

		utils::Coordinate<3> pos{1, 1, 1};

		// initialize field: no ghost cells here
		BcField bcfields(properties.size);
		decltype(bcfields.size()) zero = 0;
		allscale::api::user::algorithm::pfor(zero,bcfields.size(),[&](const auto& pos){
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
		allscale::api::user::algorithm::pfor(zero,fields.size(),[&](const auto& pos){
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
		properties.size = { 2, 2, 2 };

	    // use center position
		utils::Coordinate<3> pos{1, 1, 1};

		// initialize field
		BcField bcfields(properties.size);
		decltype(bcfields.size()) zero = 0;
		allscale::api::user::algorithm::pfor(zero,bcfields.size(),[&](const auto& pos){
			bcfields[pos].Bc = { 1.0, 1.0, 1.0 };
		});

		Field fields(properties.size + coordinate_type(1));
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
		allscale::api::user::algorithm::pfor(zero,fields.size(),[&](const auto& pos){
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

		// this test verifies E and B fields
		std::string path = std::string(PATH_TO_INPUTS) + "/tiny.inp";
		auto params = Parameters(path);

		// initialize initial properties
		InitProperties initProperties = InitProperties(params);

		// initialize universe properties
		UniverseProperties universeProperties = UniverseProperties(params);

		// initialize field
		Field fields = initFields(initProperties, universeProperties);
		auto start = utils::Coordinate<3>(1);
		auto end = universeProperties.size + coordinate_type(2);

	    // apply static field solver and check results
		allscale::api::user::algorithm::pfor(start, end, [&](const utils::Coordinate<3>& pos) {
			solveFieldStatically(universeProperties, pos, fields);
		});

		// verify for the Dipole case
		utils::Coordinate<3> pos{1, 2, 4};

		EXPECT_NEAR(fields[pos].Bext.x, -2.74633e-08, 1e-10);
		EXPECT_NEAR(fields[pos].Bext.y, -2.40304e-08, 1e-10);
		EXPECT_NEAR(fields[pos].Bext.z, 1.44183e-08, 1e-10);
		EXPECT_NEAR(fields[pos].B.x, 0.0, 1e-10);
		EXPECT_NEAR(fields[pos].B.y, 0.0, 1e-10);
		EXPECT_NEAR(fields[pos].B.z, 0.0001, 1e-10);

		allscale::api::user::algorithm::pfor(start, end, [&](const utils::Coordinate<3>& cur) {
			EXPECT_NEAR(fields[cur].E.x, -0.0, 1e-06);
			EXPECT_NEAR(fields[cur].E.y, 2.0e-06, 1e-08);
			EXPECT_NEAR(fields[cur].E.z, -0.0, 1e-06);
		});
    }


	TEST(Field, updateFieldsOnBoundaries) {
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 10, 10, 10 };
		properties.cellWidth = { 1, 1, 1 };
		properties.useCase = UseCase::Dipole;

		utils::Coordinate<3> pos{0, 0, 0};

		// initialize field
		Field field(properties.size + coordinate_type(3));
		decltype(field.size()) start = 1;
		allscale::api::user::algorithm::pfor(start, field.size() - coordinate_type(1), [&](const auto& pos){
			field[pos].E = { 1.0, 1.0, 1.0 };
			field[pos].B = { 2.0, 2.0, 2.0 };
		});
		BcField bcfield(properties.size + coordinate_type(2));
		allscale::api::user::algorithm::pfor(start, bcfield.size() - coordinate_type(1), [&](const auto& pos){
			bcfield[pos].Bc = { 3.0, 3.0, 3.0 };
		});

		updateFieldsOnBoundaries(field, bcfield);

		auto E = field[ {0, 5, 11} ].E;
		EXPECT_NEAR( E.x, 1.0, 1e-06 );
		EXPECT_NEAR( E.y, 1.0, 1e-06 );
		EXPECT_NEAR( E.z, 1.0, 1e-06 );

		auto B = field[ {9, 12, 2} ].B;		
		EXPECT_NEAR( B.x, 2.0, 1e-06 );
		EXPECT_NEAR( B.y, 2.0, 1e-06 );
		EXPECT_NEAR( B.z, 2.0, 1e-06 );

		auto Bc = bcfield[ {3, 8, 12} ].Bc;
		EXPECT_NEAR( Bc.x, 3.0, 1e-06 );
		EXPECT_NEAR( Bc.y, 3.0, 1e-06 );
		EXPECT_NEAR( Bc.z, 3.0, 1e-06 );

		allscale::api::user::algorithm::pfor(start, field.size() - coordinate_type(1), [&](const auto& pos){
			field[pos].E = { (double) pos[0], (double) pos[1], (double) pos[2] };
			field[pos].B = { 2.0 * pos[0], 2.0 * pos[1], 2.0 * pos[2]};
		});

		updateFieldsOnBoundaries(field, bcfield);

		E = field[ {0, 5, 11} ].E;
		EXPECT_NEAR( E.x, 11.0, 1e-06 );
		EXPECT_NEAR( E.y, 5.0, 1e-06 );
		EXPECT_NEAR( E.z, 11.0, 1e-06 );

		B = field[ {9, 12, 2} ].B;		
		EXPECT_NEAR( B.x, 18.0, 1e-06 );
		EXPECT_NEAR( B.y, 2.0, 1e-06 );
		EXPECT_NEAR( B.z, 4.0, 1e-06 );
	}


	TEST(Field, solveFieldForward) {
		
		// this test verifies the explicit forward method§§§
		// for a moment the test as well as the code are empty
		
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 5,5,5 };
		properties.cellWidth = { 1.0,1.0,1.0 };
		properties.dt = 1.0;

		utils::Coordinate<3> pos{2, 2, 2};
		properties.useCase = UseCase::Dipole;

		// initialize field
		Field field(properties.size + coordinate_type(3));
		decltype(field.size()) zero = 0;
		allscale::api::user::algorithm::pfor(zero,field.size(),[&](const auto& pos){
			field[pos].E = { 0.2, 0.2, 0.2 };
			field[pos].B = { 0.4, 0.4, 0.4 };
		});
		// initialize field
		BcField bcfield(properties.size + coordinate_type(2));
		allscale::api::user::algorithm::pfor(zero,bcfield.size(),[&](const auto& pos){
			bcfield[pos].Bc = { 0.8, 0.8, 0.8 };
		});
		// initialize density
		CurrentDensity density(properties.size + coordinate_type(1));
		allscale::api::user::algorithm::pfor(zero,density.size(),[&](const auto& pos){
			density[pos].J = { 0.6, 0.6, 0.6 };
		});
		

	    // apply leapfrog field solver and check results
		solveFieldForward(properties, pos, density, field, bcfield);
		auto E = field[pos].E;
		EXPECT_NEAR( E.x, -0.4, 1e-15 );
		EXPECT_NEAR( E.y, -0.4, 1e-15 );
		EXPECT_NEAR( E.z, -0.4, 1e-15 );

		auto Bc = bcfield[pos].Bc;
		EXPECT_NEAR( Bc.x, 0.8, 1e-15 );
		EXPECT_NEAR( Bc.y, 0.8, 1e-15 );
		EXPECT_NEAR( Bc.z, 0.8, 1e-15 );
	}


	TEST(Field, solveFieldForwardTrajectory) {
		
		// this test verifies the explicit forward method
		// for a moment the test as well as the code are empty
		
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 5,5,5 };
		properties.cellWidth = { 1.0,1.0,1.0 };
		properties.dt = 0.1;

		properties.useCase = UseCase::Dipole;

		// initialize field
		Field field(properties.size + coordinate_type(3));
		decltype(field.size()) zero = 0;
		allscale::api::user::algorithm::pfor(zero,field.size(),[&](const auto& pos){
			field[pos].E = { 0.2, 0.2, 0.2 };
			field[pos].B = { 0.4, 0.4, 0.4 };
		});
		// initialize bcfield
		BcField bcfield(properties.size + coordinate_type(2));
		allscale::api::user::algorithm::pfor(zero,bcfield.size(),[&](const auto& pos){
			bcfield[pos].Bc = { 0.8, 0.8, 0.8 };
		});
		// initialize density
		CurrentDensity density(properties.size + coordinate_type(1));
		allscale::api::user::algorithm::pfor(zero,density.size(),[&](const auto& pos){
			density[pos].J = { 0.6, 0.6, 0.6 };
		});
		

	    // apply the field solver and check results
		unsigned numSteps = 20;
		utils::Coordinate<3> field_pos{2, 2, 2};
		for(unsigned i = 0; i < numSteps; ++i) {
			solveFieldForward(properties, field_pos, density, field, bcfield);
			auto E = field[field_pos].E;
			auto B = field[field_pos].B;
			auto Bc = bcfield[field_pos].Bc;

			std::cout << i << " " << i * properties.dt << " ";
			std::cout << E.x << " " << E.y << " " << E.z << " ";
			std::cout << Bc.x << " " << Bc.y << " " << Bc.z << " ";
			std::cout << B.x << " " << B.y << " " << B.z << "\n";
		}
	}


	TEST(Field, solveFieldLeapfrog) {
		
		// this test verifies the leapfrog algorithm
		// for a moment the test as well as the code are empty
		
		// Set universe properties
		UniverseProperties properties;
		properties.size = { 5,5,5 };
		properties.cellWidth = { 1,1,1 };
		properties.useCase = UseCase::Dipole;

		utils::Coordinate<3> pos{2, 2, 2};

		// initialize field
		Field field(properties.size + coordinate_type(1));
		decltype(field.size()) zero = 0;
		allscale::api::user::algorithm::pfor(zero,field.size(),[&](const auto& pos){
			field[pos].E = { 0.0, 0.0, 0.0 };
			field[pos].B = { 0.0, 0.0, 0.0 };
		});
		// initialize bcfield
		BcField bcfield(properties.size + coordinate_type(2));
		allscale::api::user::algorithm::pfor(zero,bcfield.size(),[&](const auto& pos){
			bcfield[pos].Bc = { 0.0, 0.0, 0.0 };
		});
		// initialize density
		CurrentDensity density(properties.size + coordinate_type(1));
		allscale::api::user::algorithm::pfor(zero,density.size(),[&](const auto& pos){
			density[pos].J = { 0.0, 0.0, 0.0 };
		});

	    // apply leapfrog field solver and check results
		solveFieldLeapfrog(properties, pos, density, field, bcfield);
		auto E = field[pos].E;
		EXPECT_NEAR( E.x, 0.0, 1e-15 );
		EXPECT_NEAR( E.y, 0.0, 1e-15 );
		EXPECT_NEAR( E.z, 0.0, 1e-15 );

		auto B = field[pos].B;
		EXPECT_NEAR( B.x, 0.0, 1e-15 );
		EXPECT_NEAR( B.y, 0.0, 1e-15 );
		EXPECT_NEAR( B.z, 0.0, 1e-15 );
	}

} // end namespace ipic3d

