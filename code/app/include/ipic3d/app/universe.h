#pragma once

#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"

#include "allscale/api/user/operator/pfor.h"

namespace ipic3d {

	/**
	* A structure representing the simulated Universe including both cells and the field of forces
	*/
	struct Universe {

		// The cells of this universe
		Cells cells;

		// The field of this universe
		Field field;

		using coordinate_type = allscale::api::user::data::GridPoint<3>;

		/**
		* Creates a Universe of cells and a field of forces the given size.
		* The dimensions of the field grid exceed the dimensions of the cell grid by 1 in every dimension.
		*
		* @param dims the size of the Universe (equal to the size of the grid of cells)
		*/
		Universe(const coordinate_type& dims) : cells(Cells(dims)), field(Field(dims+coordinate_type(1))) {}

		Universe(Cells&& cells, Field&& field) : cells(std::move(cells)), field(std::move(field)) {}

		coordinate_type size() const {
			return cells.size();
		}

	};

}
