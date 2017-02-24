#pragma once

#include "allscale/api/user/data/grid.h"

#include "ipic3d/app/vector.h"

namespace ipic3d {

	using coordinate_type = allscale::api::user::data::GridPoint<3>;

	/**
	 * Holds all properties common to/uniform in a universe
	 *
	*/
	struct UniverseProperties {

		// The size of this universe
		coordinate_type size;

		// The width of cells
		Vector3<double> cellWidth;

		UniverseProperties(const coordinate_type& size = { 1,1,1 }, const Vector3<double>& cellWidth = {1,1,1}) : size(size), cellWidth(cellWidth) {
			assert_true(size.x > 0 && size.y > 0 && size.z > 0) << "Expected positive non-zero universe size, but got " << size;
			assert_true(cellWidth.x > 0 && cellWidth.y > 0 && cellWidth.z > 0) << "Expected positive non-zero cell widths, but got " << cellWidth;
		}

		friend std::ostream& operator<<(std::ostream& out, const UniverseProperties& props) {
			out << "Universe properties:" << std::endl;
			out << "\tSize: " << props.size << std::endl;
			out << "\tCell width: " << props.cellWidth << std::endl;
			return out;
		}

	};

	Vector3<double> getCenterOfCell(const coordinate_type& pos, const UniverseProperties& properties) {
		Vector3<double> tempPos{ (double)pos.x, (double)pos.y, (double)pos.z };
		return elementwiseProduct(tempPos, properties.cellWidth) + properties.cellWidth / 2;
	}

}
