#pragma once

#include "ipic3d/app/cell.h"
#include "ipic3d/app/field.h"
#include "ipic3d/app/universe_properties.h"

namespace ipic3d {

	/**
	* A structure representing the simulated Universe including both cells and the field of forces
	*/
	struct Universe {

		// Uniform properties of this universe
		const UniverseProperties properties;

		// The cells of this universe
		Cells cells;

		// The field of this universe
		Field field;

		// The magnetic field of this universe defined on centers of cells
		BcField bcfield;

		/**
		* Creates a Universe of cells and fields of forces of the given size.
		* The dimensions of the field grid exceed the dimensions of the cell grid by 1 in every dimension.
		*
		* @param properties the properties of the universe to be created
		*/
	    Universe(const UniverseProperties& properties = UniverseProperties())
	        : properties(properties), cells(Cells(properties.size)), field(Field(properties.size + coordinate_type(3))), bcfield(BcField(properties.size + coordinate_type(2))) // two for the two extra boundary cells and one as fields are defined on nodes of the cells
		{ 
			auto dims = properties.size;
			assert_true(dims.x > 0 || dims.y > 0 || dims.z > 0) << "Expected positive non-zero dimensions, but got " << dims;
		}

	    Universe(const UniverseProperties& properties, Cells&& cells, Field&& field, BcField&& bcfield) : properties(properties), cells(std::move(cells)), field(std::move(field)), bcfield(std::move(bcfield)) {
			auto size = cells.size();
			assert_true(size == properties.size) << "Expected size of universe and size of cell grid to match, but got " << size << " and " << properties.size;
			assert_true((size + coordinate_type(3)) == field.size()) << "Expected size of field grid to be equal to size of cell grid + 1 but got " << field.size() << " and " << size;
			assert_true(size.x > 0 || size.y > 0 || size.z > 0) << "Expected positive non-zero dimensions, but got " << size;
		}

	};

	Universe createUniverseFromParams(const Parameters& params) {

		// initialize initial properties
		InitProperties initProperties = InitProperties(params);

		std::cout << initProperties;

		// initialize universe properties
		UniverseProperties universeProperties = UniverseProperties(params);

		std::cout << universeProperties;

		// initialize grid of cells
		Cells&& cells = initCells(initProperties, universeProperties);

		// initialize fields on node
		Field&& field = initFields(initProperties, universeProperties);

		// initialize magnetic fields on centers
		BcField&& bcField = initBcFields(universeProperties, field);

		// create the universe with the given properties, cells and fields
		Universe universe(universeProperties, std::move(cells), std::move(field), std::move(bcField));

		return universe;
	}

}
