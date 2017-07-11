#include <gtest/gtest.h>

#include "ipic3d/app/cell.h"
#include "ipic3d/app/universe.h"

namespace ipic3d {

	TEST(Cell, initCells) {

		// this test checks the initCells() function

		std::string path = std::string(PATH_TO_INPUTS) + "/micro.inp";
		auto params = Parameters::read(path);

		// initialize initial properties
		InitProperties initProperties = InitProperties(params);

//		Cells cells = initCells(initProperties, properties);
		
		Universe universe = createUniverseFromParams(params);	

		// verify the number of particles per cell
		int particlesPerCell = initProperties.particlesPerCell[0].x * initProperties.particlesPerCell[0].y * initProperties.particlesPerCell[0].z;

		utils::Coordinate<3> zero = 0;
		allscale::api::user::pfor(zero, universe.properties.size, [&](const utils::Coordinate<3>& pos) {
	    	EXPECT_EQ(particlesPerCell, (int) universe.cells[pos].particles.size());		
		});
	}


	TEST(Cell, VerifyCorrectParticlesPositionInCell) {
	
		// this test checkes the correct placement of particles within cells

		// Set universe properties
		UniverseProperties properties;
		properties.size = {4,4,4};
		properties.cellWidth = { .5,.5,.5 };
		properties.dt = 0.1;
		properties.useCase = UseCase::Test;

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
		p0.position.x = 0.4; p0.position.y = p0.position.z = 0.4;
		p0.velocity.x = 0.6; p0.velocity.y = p0.velocity.z = 0.6;
		p0.q = p0.mass = 1;
		p1.position.x = 0.8; p1.position.y = p1.position.z = 0.4;
		p1.velocity.x = 0.3; p1.velocity.y = p1.velocity.z = 0.3;
		p1.q = p1.mass = 1;
		p2.position.y = 0.8; p2.position.x = p2.position.z = 0.4;
		p2.velocity.x = 0.1; p2.velocity.y = p2.velocity.z = 0.1;
		p2.q = p2.mass = 1;
		p3.position.x = 0.75; p3.position.y = 0.75; p3.position.z = 0.4;
		p3.velocity.x = 0.8; p3.velocity.y = p3.velocity.z = 0.8;
		p3.q = p3.mass = 1;
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
		ASSERT_FALSE( VerifyCorrectParticlesPositionInCell(properties, a, {0,0,0}) );
		ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, b, {1,0,0}) );
		ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, c, {0,1,0}) );
		ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, d, {1,1,0}) );
		ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, e, {0,0,1}) );
		ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, f, {1,0,1}) );
		ASSERT_TRUE( VerifyCorrectParticlesPositionInCell(properties, g, {0,1,1}) );
		ASSERT_FALSE( VerifyCorrectParticlesPositionInCell(properties, h, {1,1,1}) );
		
	}

}
