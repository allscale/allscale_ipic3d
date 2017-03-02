#pragma once

#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/universe_properties.h"

namespace ipic3d {

	/**
	* A structure representing the simulated Universe including both cells and the field of forces
	*/
	struct Universe {

		// The cells of this universe
		Cells cells;

		// The field of this universe
		Field field;

		// Uniform properties of this universe
		const UniverseProperties properties;

		/**
		* Creates a Universe of cells and a field of forces the given size.
		* The dimensions of the field grid exceed the dimensions of the cell grid by 1 in every dimension.
		*
		* @param dims the size of the Universe (equal to the size of the grid of cells)
		*/
	    Universe(const UniverseProperties& properties = UniverseProperties())
	        : cells(Cells(properties.size)), field(Field(properties.size + coordinate_type(1))), properties(properties) {
			auto dims = properties.size;
			assert_true(dims.x > 0 || dims.y > 0 || dims.z > 0) << "Expected positive non-zero dimensions, but got " << dims;
		}

	    Universe(Cells&& cells, Field&& field, const UniverseProperties& properties) : cells(std::move(cells)), field(std::move(field)), properties(properties) {
			auto dims = cells.size();
			assert_true(dims == properties.size) << "Expected size of universe and size of cell grid to match, but got " << dims << " and " << properties.size;
			assert_true(dims.x > 0 || dims.y > 0 || dims.z > 0) << "Expected positive non-zero dimensions, but got " << dims;
		}

	};

}
