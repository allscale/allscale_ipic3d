#pragma once

#include "ipic3d/app/universe.h"

namespace ipic3d{

	int countParticlesInDomain(Universe& universe) {

		// get the number of particles in all cells before the simulation begins
		std::atomic<int> particles = ATOMIC_VAR_INIT(0);;
		decltype(universe.field.size()) zero = 0;
		allscale::api::user::pfor(zero,universe.properties.size,[&](auto& pos){
			std::atomic_fetch_add( &particles, (int) universe.cells[pos].particles.size() );
		});

		return particles;
	}

}
