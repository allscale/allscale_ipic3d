#pragma once

#include <array>
#include <initializer_list>

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/data/vector.h"

namespace ipic3d {
namespace utils {

	template<typename T, size_t Dims>
	using vec = allscale::api::user::data::Vector<T,Dims>;

	template<size_t Dims>
	using Coordinate = allscale::api::user::data::GridPoint<Dims>;

	template<size_t Dims>
	using Size = Coordinate<Dims>;


} // end namespace utils
} // end namespace ipic3d
