#include <gtest/gtest.h>

#include "ipic3d/app/cell.h"
#include "ipic3d/app/universe.h"

#include "allscale/api/core/io.h"

namespace ipic3d {

	TEST(Cell, initCells) {

		// this test checks the initCells() function

		std::string path = std::string(PATH_TO_INPUTS) + "/test.inp";
		auto params = Parameters(path);

		// initialize initial properties
		InitProperties initProperties = InitProperties(params);
		UniverseProperties universeProperties = UniverseProperties(params);

		Cells&& cells = initCells(params, initProperties, universeProperties);

		// verify the number of particles per cell
		int particlesPerCell = initProperties.particlesPerCell[0].x * initProperties.particlesPerCell[0].y * initProperties.particlesPerCell[0].z;

		utils::Coordinate<3> zero = 0;
		allscale::api::user::algorithm::pfor(zero, universeProperties.size, [&](const utils::Coordinate<3>& pos) {
			EXPECT_EQ(particlesPerCell, (int) cells[pos].particles.size());
		});
	}


	TEST(Cell, VerifyCorrectParticlesPositionInCell) {
	
		// this test checks the correct placement of particles within cells

		// Set universe properties
		UniverseProperties properties;
		properties.size = {4,4,4};
		properties.cellWidth = { .5,.5,.5 };
		properties.dt = 0.1;
		properties.useCase = UseCase::Dipole;

		// Create a universe with these properties
		Universe universe = Universe(properties);

		// configure the cell
		Cell& a = universe.cells[{0,0,0}];
		Cell& b = universe.cells[{1,0,0}];
		Cell& c = universe.cells[{0,1,0}];
		Cell& d = universe.cells[{1,1,0}];
		Cell& e = universe.cells[{0,0,1}];
		Cell& f = universe.cells[{1,0,1}];
		Cell& g = universe.cells[{0,1,1}];
		Cell& h = universe.cells[{1,1,1}];

		// add particles: one per cells in the subdomain {2,2,2}
		Particle p0, p1, p2, p3, p4, p5, p6, p7;
		p0.position.x = p0.position.y = p0.position.z = 0.4;
		p0.velocity.x = p0.velocity.y = p0.velocity.z = 0.6;
		p0.q = p0.qom = 1.0;
		p1.position.x = 0.8; p1.position.y = p1.position.z = 0.4;
		p1.velocity.x = p1.velocity.y = p1.velocity.z = 0.3;
		p1.q = p1.qom = 1.0;
		p2.position.y = 0.8; p2.position.x = p2.position.z = 0.4;
		p2.velocity.x = p2.velocity.y = p2.velocity.z = 0.1;
		p2.q = p2.qom = 1.0;
		p3.position.x = p3.position.y = 0.75; p3.position.z = 0.4;
		p3.velocity.x = p3.velocity.y = p3.velocity.z = 0.8;
		p3.q = p3.qom = 1.0;
		p4 = p0;
		p4.position.z = 0.6;
		p5 = p1;
		p5.position.z = 0.8;
		p6 = p3;
		p6.position.x = 0.25;
		p6.position.z = 0.75;
		p7 = p3;
		p7.position.z = 0.8;

		// add incorrectly placed particles in two cells
		Particle p01, p02, p03, p71, p72;
		p01 = p0;
		p01.position.x = 0.6;
		p02 = p0;
		p02.position.y = 0.6;
		p03 = p0;
		p03.position.z = 0.6;
		p71 = p7;
		p71.position.x = 1.75;
		p72 = p7;
		p72.position.y = p72.position.z = 1.08;

		// add test particle to first cell
		a.particles.push_back(p0);
		a.particles.push_back(p01);
		a.particles.push_back(p02);
		a.particles.push_back(p03);
		b.particles.push_back(p1);
		c.particles.push_back(p2);
		d.particles.push_back(p3);
		e.particles.push_back(p4);
		f.particles.push_back(p5);
		g.particles.push_back(p6);
		h.particles.push_back(p7);
		h.particles.push_back(p7);
		h.particles.push_back(p71);
		h.particles.push_back(p72);

		// verify the proper placement of particles 
		ASSERT_FALSE( verifyCorrectParticlesPositionInCell(properties, a, {0,0,0}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, b, {1,0,0}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, c, {0,1,0}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, d, {1,1,0}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, e, {0,0,1}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, f, {1,0,1}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, g, {0,1,1}) );
		ASSERT_FALSE( verifyCorrectParticlesPositionInCell(properties, h, {1,1,1}) );
	}


	TEST(Cell, projectToDensityField) {
	
		// this test verifies projection to density fields defined on nodes

		// Set universe properties
		UniverseProperties properties;
		properties.size = {2,2,2};
		properties.cellWidth = { .5,.5,.5 };
		properties.dt = 0.1;
		properties.useCase = UseCase::Dipole;

		// Create a universe with these properties
		Universe universe = Universe(properties);

		auto zero = utils::Coordinate<3>(0);

		// init current density
		allscale::api::user::algorithm::pfor(zero, universe.currentDensity.size(), [&](const utils::Coordinate<3>& pos) {
			universe.currentDensity[pos].J = {0.0,0.0,0.0};
		});

		// configure the cell
		Cell& a = universe.cells[{0,0,0}];
		Cell& b = universe.cells[{1,0,0}];
		Cell& c = universe.cells[{0,1,0}];
		Cell& d = universe.cells[{1,1,0}];
		Cell& e = universe.cells[{0,0,1}];
		Cell& f = universe.cells[{1,0,1}];
		Cell& g = universe.cells[{0,1,1}];
		Cell& h = universe.cells[{1,1,1}];

		// add particles: one per cells in the subdomain {2,2,2}
		Particle p0, p1, p2, p3, p4, p5, p6, p7;
		p0.position.x = p0.position.y = p0.position.z = 0.4;
		p0.velocity.x = p0.velocity.y = p0.velocity.z = 0.6;
		p0.q = p0.qom = 1.0;
		p1.position.x = 0.8; p1.position.y = p1.position.z = 0.4;
		p1.velocity.x = p1.velocity.y = p1.velocity.z = 0.3;
		p1.q = p1.qom = 1.0;
		p2.position.y = 0.8; p2.position.x = p2.position.z = 0.4;
		p2.velocity.x = p2.velocity.y = p2.velocity.z = 0.1;
		p2.q = p2.qom = 1.0;
		p3.position.x = p3.position.y = 0.75; p3.position.z = 0.4;
		p3.velocity.x = p3.velocity.y = p3.velocity.z = 0.8;
		p3.q = p3.qom = 1.0;
		p4 = p0;
		p4.position.z = 0.6;
		p5 = p1;
		p5.position.z = 0.8;
		p6 = p3;
		p6.position.x = 0.25;
		p6.position.z = 0.75;
		p7 = p3;
		p7.position.z = 0.8;

		// add test particle to first cell
		a.particles.push_back(p0);
		b.particles.push_back(p1);
		c.particles.push_back(p2);
		d.particles.push_back(p3);
		e.particles.push_back(p4);
		f.particles.push_back(p5);
		g.particles.push_back(p6);
		h.particles.push_back(p7);

		// verify the proper placement of particles 
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, a, {0,0,0}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, b, {1,0,0}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, c, {0,1,0}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, d, {1,1,0}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, e, {0,0,1}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, f, {1,0,1}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, g, {0,1,1}) );
		ASSERT_TRUE( verifyCorrectParticlesPositionInCell(properties, h, {1,1,1}) );
	
		auto size = universe.cells.size();

		allscale::api::user::data::Grid<DensityNode, 3> densityContributions(size * 2);

		// compute current density
		allscale::api::user::algorithm::pfor(zero, size, [&](const utils::Coordinate<3>& pos) {
			projectToDensityField(properties, universe.cells[pos], pos, densityContributions);
		});

		allscale::api::user::algorithm::pfor(zero, size, [&](const utils::Coordinate<3>& pos) {
			aggregateDensityContributions(properties, densityContributions, pos, universe.currentDensity[pos]);
		});

		// verify computed results
		auto Js = universe.currentDensity[utils::Coordinate<3>({1,1,1})].J;
		EXPECT_NEAR( Js.x, 1.0952, 1e-4);	
		EXPECT_NEAR( Js.y, 1.0952, 1e-4);	
		EXPECT_NEAR( Js.z, 1.0952, 1e-4);	
	}
		
	
	TEST(Cell, TestCellOutput) {

		// this test checks the output of the number of particles per cell

		// Set universe properties
		UniverseProperties properties;
		properties.size = { 2,2,2 };
		properties.cellWidth = { .5,.5,.5 };
		properties.dt = 0.1;
		properties.useCase = UseCase::Dipole;

		// Create a universe with these properties
		Universe universe = Universe(properties);

		// configure the cell
		Cell& b = universe.cells[{1, 0, 0}];
		Cell& c = universe.cells[{0, 1, 0}];
		Cell& d = universe.cells[{1, 1, 0}];
		Cell& e = universe.cells[{0, 0, 1}];
		Cell& f = universe.cells[{1, 0, 1}];
		Cell& g = universe.cells[{0, 1, 1}];
		Cell& h = universe.cells[{1, 1, 1}];

		// create dummy particle
		Particle p;

		// add dummy particle to cells for counting only
		b.particles.push_back(p);
		c.particles.push_back(p);
		c.particles.push_back(p);
		d.particles.push_back(p);
		d.particles.push_back(p);
		d.particles.push_back(p);
		e.particles.push_back(p);
		f.particles.push_back(p);
		g.particles.push_back(p);
		h.particles.push_back(p);
		h.particles.push_back(p);
		h.particles.push_back(p);
		h.particles.push_back(p);

		// verify the proper particles count output
		allscale::api::core::BufferIOManager manager;
		auto text = manager.createEntry("text");
		auto out = manager.openOutputStream(text);
		outputNumberOfParticlesPerCell(universe.cells, out);
		manager.close(out);

		auto in = manager.openInputStream(text);
		coordinate_type size;
		char separator;

		// read size
		EXPECT_TRUE(in >> separator);
		EXPECT_TRUE(in >> size.x);
		EXPECT_TRUE(in >> separator);
		EXPECT_TRUE(in >> size.y);
		EXPECT_TRUE(in >> separator);
		EXPECT_TRUE(in >> size.z);
		EXPECT_TRUE(in >> separator);

		EXPECT_EQ((coordinate_type{ 2,2,2 }), size);

		// read payload
		for(int i = 0; i < size.x * size.y * size.z; ++i) {
			// read index
			coordinate_type index;
			EXPECT_TRUE(in >> index.x) << "Trying to read packet " << i;
			EXPECT_TRUE(in >> separator) << "Trying to read packet " << i;
			EXPECT_TRUE(in >> index.y) << "Trying to read packet " << i;
			EXPECT_TRUE(in >> separator) << "Trying to read packet " << i;
			EXPECT_TRUE(in >> index.z) << "Trying to read packet " << i;
			EXPECT_TRUE(in >> separator) << "Trying to read packet " << i;

			// read and verify data
			std::size_t particleCount;
			EXPECT_TRUE(in >> particleCount) << "Trying to read packet " << i;
			EXPECT_EQ(particleCount, universe.cells[index].particles.size());
		}

	}

	TEST(Cell, DISABLED_DensityContributions) {
		// TODO: implement test for density contribution computation
	}

}
