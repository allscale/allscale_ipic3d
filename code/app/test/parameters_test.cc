#include <gtest/gtest.h>

#include "ipic3d/app/parameters.h"

namespace ipic3d {
	
	TEST(Parameters, Parser) {

		// this test verifies the parser

		// test wrong file name
		std::string path = std::string(PATH_TO_INPUTS) + "/test.inp1";
		auto params = Parameters(path);

		// test with the correct path
		path = std::string(PATH_TO_INPUTS) + "/test.inp";
		params = Parameters(path);

		EXPECT_NEAR(params.c, 1.0, 1e-15);
		EXPECT_NEAR(params.dt, 0.15, 1e-15);
		EXPECT_EQ(params.ncycles, 10);

		EXPECT_NEAR(params.L.x, 10.0, 1e-15);
		EXPECT_NEAR(params.L.y, 10.0, 1e-15);
		EXPECT_NEAR(params.L.z, 10.0, 1e-15);

		EXPECT_NEAR(params.objectCenter.x, 6.0, 1e-15);
		EXPECT_NEAR(params.objectCenter.y, 5.0, 1e-15);
		EXPECT_NEAR(params.objectCenter.z, 5.0, 1e-15);
		EXPECT_NEAR(params.planetRadius, 0.5, 1e-15);
		EXPECT_NEAR(params.delta, 0.5, 1e-15);

		EXPECT_EQ(params.ncells.x, 4);
		EXPECT_EQ(params.ncells.y, 4);
		EXPECT_EQ(params.ncells.z, 4);

		EXPECT_NEAR(params.dspace.x, 2.5, 1e-15);
		EXPECT_NEAR(params.dspace.y, 2.5, 1e-15);
		EXPECT_NEAR(params.dspace.z, 2.5, 1e-15);

		EXPECT_EQ(params.ns, 2);
		for (auto it = params.npcelx.cbegin(); it != params.npcelx.cend(); ++it) 
			EXPECT_EQ(*it, 5);
		for (auto it = params.npcely.cbegin(); it != params.npcely.cend(); ++it) 
			EXPECT_EQ(*it, 5);
		for (auto it = params.npcelz.cbegin(); it != params.npcelz.cend(); ++it) 
			EXPECT_EQ(*it, 5);

		EXPECT_NEAR(params.qom[0], -25.0, 1e-15);
		EXPECT_NEAR(params.qom[1], 1.0, 1e-15);

		EXPECT_NEAR(params.rhoINIT[0], 1.0, 1e-15);
		EXPECT_NEAR(params.rhoINIT[1], 1.0, 1e-15);

		EXPECT_NEAR(params.uth[0], 0.045, 1e-15);
		EXPECT_NEAR(params.uth[1], 0.0063, 1e-15);

		EXPECT_NEAR(params.vth[0], 0.045, 1e-15);
		EXPECT_NEAR(params.vth[1], 0.0063, 1e-15);

		EXPECT_NEAR(params.wth[0], 0.045, 1e-15);
		EXPECT_NEAR(params.wth[1], 0.0063, 1e-15);

		EXPECT_NEAR(params.u0[0], 0.02, 1e-15);
		EXPECT_NEAR(params.u0[1], 0.02, 1e-15);

		EXPECT_NEAR(params.v0[0], 0.0, 1e-15);
		EXPECT_NEAR(params.v0[1], 0.0, 1e-15);

		EXPECT_NEAR(params.w0[0], 0.0, 1e-15);
		EXPECT_NEAR(params.w0[1], 0.0, 1e-15);

		EXPECT_TRUE( params.useCase == UseCase::Dipole );
	
		EXPECT_TRUE( params.wmethod.compare("pvtk") == 0 );
		EXPECT_TRUE( params.SimName.compare("Dipole3D") == 0 );
		EXPECT_TRUE( params.SaveDirName.compare("data") == 0 );

		EXPECT_TRUE( params.PoissonCorrection.compare("no") == 0 );

		EXPECT_NEAR(params.B0.x, 0.0, 1e-15);
		EXPECT_NEAR(params.B0.y, 0.0, 1e-15);
		EXPECT_NEAR(params.B0.z, 0.0001, 1e-15);

		EXPECT_NEAR(params.B1.x, 0.0, 1e-15);
		EXPECT_NEAR(params.B1.y, 0.0, 1e-15);
		EXPECT_NEAR(params.B1.z, 2.0, 1e-15);

		EXPECT_EQ(params.FieldOutputCycle, 1);
		EXPECT_TRUE( params.FieldOutputTag.compare("B+E+Je+Ji") == 0 );
		EXPECT_TRUE( params.MomentsOutputTag.compare("rho+PXX+PXY+PXZ+PYY+PYZ+PZZ") == 0 );

		EXPECT_EQ(params.ParticlesOutputCycle, 10);
		EXPECT_TRUE( params.ParticlesOutputTag.compare("position+velocity+q") == 0 );
	}

} // end namespace ipic3d
