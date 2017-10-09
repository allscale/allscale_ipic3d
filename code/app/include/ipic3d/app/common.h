#pragma once

#include "ipic3d/app/universe.h"

#include "allscale/api/user/operator/preduce.h"

namespace ipic3d{

	int countParticlesInDomain(Universe& universe) {

		// get the number of particles in all cells before the simulation begins
		
		// TODO: use more convenient reduction operators once they are available in the API
		auto map = [&](const coordinate_type& index, int& res) {
			res += (int)universe.cells[index].particles.size();
		};

		auto reduce = [&](const int& a, const int& b) { return a + b; };
		auto init = []() { return 0; };

		coordinate_type zero(0);
		coordinate_type full(universe.cells.size());

		return allscale::api::user::preduce(zero, full, map, reduce, init).get();

	}

}
