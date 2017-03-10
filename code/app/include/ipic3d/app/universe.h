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

		// The magnetic field of this universe defined on centers of cells
		BcField bcfield;	

		// Uniform properties of this universe
		UniverseProperties properties;

		/**
		* Creates a Universe of cells and a field of forces the given size.
		* The dimensions of the field grid exceed the dimensions of the cell grid by 1 in every dimension.
		*
		* @param dims the size of the Universe (equal to the size of the grid of cells)
		*/
	    Universe(const UniverseProperties& properties = UniverseProperties())
	        : cells(Cells(properties.size)), field(Field(properties.size + coordinate_type(1))), bcfield(properties.size), properties(properties) {
			auto dims = properties.size;
			assert_true(dims.x > 0 || dims.y > 0 || dims.z > 0) << "Expected positive non-zero dimensions, but got " << dims;
		}

	    Universe(Cells&& cells, Field&& field, BcField&& bcfield, const UniverseProperties& properties) : cells(std::move(cells)), field(std::move(field)), bcfield(std::move(bcfield)), properties(properties) {
			auto size = cells.size();
			assert_true(size == properties.size) << "Expected size of universe and size of cell grid to match, but got " << size << " and " << properties.size;
			assert_true((size + coordinate_type(1)) == field.size()) << "Expected size of field grid to be equal to size of cell grid + 1 but got " << field.size() << " and " << size;
			assert_true(size.x > 0 || size.y > 0 || size.z > 0) << "Expected positive non-zero dimensions, but got " << size;
		}

		Universe(const Parameters& params) {

			// initialize initial properties
			InitProperties initProperties = InitProperties(params);

			std::cout << initProperties;

			// initialize universe properties
			properties = UniverseProperties(params);

			std::cout << properties;

			//Universe universe(initCells(initProperties, universeProperties), initFields(initProperties, universeProperties), initBcFields(universeProperties), universeProperties);

			// initialize grid of cells
			//cells = initCells( initProperties, properties );
			cells(properties.size);

			// initialize fields on node
			//field = initFields( initProperties, properties );
			field(properties.size);

			// initialize magnetic fields on centers
			//bcfield = initBcFields( properties, field );
			bvfield(properties.size);
		}

	};

}
