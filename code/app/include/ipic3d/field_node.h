#pragma once

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/data/vector.h"

namespace ipic3d {

	template<typename T> struct Vector {
		T x;
		T y;
		T z;
	};

	struct FieldNode {
		Vector<double> E;				// electric field components defined on nodes
		Vector<double> B;				// magnetic field components defined on nodes
		Vector<double> Bc;				// magnetic field components defined on central points between nodes TODO: to clarify this
		Vector<double> Bext;			// external magnetic field on nodes
	};

	using Field = allscale::api::user::data::Grid<FieldNode,3>;	// a 3D grid of field nodes

} // end namespace ipic3d
