#pragma once

#include <cmath>
#include "allscale/utils/vector.h"

namespace ipic3d {

	template<typename T>
	using Vector3 = allscale::utils::Vector<T,3>;

	template<typename T, std::size_t Dims>
	T norm(const allscale::utils::Vector<T,Dims>& vec) {
		return std::sqrt(allscale::utils::sumOfSquares(vec));
	}

} // end namespace ipic3d
