#pragma once

namespace ipic3d {

	struct Particle {
		double x,y,z;					// position (absolute - TODO: relative to cell center)
		double dx,dy,dz;				// velocity
		double q;						// charge divided over the mass of spacies
		double vxstar, vystar, vzstar;  // auxiliary parameters for the Boris mover
	};

} // end namespace ipic3d
