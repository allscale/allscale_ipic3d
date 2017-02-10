#pragma once

#include "ipic3d/cell.h"
#include "ipic3d/field_node.h"

namespace ipic3d {

	/**
	 *
	 */
	void simulateStep(double dt, Cells& cells, Field& field);

	void simulateSteps(int numSteps, double dt, Cells& cells, Field& field);

} // end namespace ipic3d
