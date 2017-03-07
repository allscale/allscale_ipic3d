#pragma once

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/data/vector.h"

#include "ipic3d/app/universe_properties.h"
#include "ipic3d/app/utils/points.h"

namespace ipic3d {

	struct FieldNode {
		// electric field components defined on nodes
		Vector3<double> E;
		// magnetic field components defined on nodes
		// TODO: should possibly be defined on its own grid with nodes inbetween the electric field nodes
		Vector3<double> B;
		// external magnetic field from dipole defined on nodes
		// TODO: read-only should only influence universe but not be modified by field solver
		Vector3<double> Bext;
	};

	// a 3D grid of field nodes
	using Field = allscale::api::user::data::Grid<FieldNode,3>;

	/**
	* TODO: provide a reference and a description of how the static field solver works
	*/
	void solveFieldStatically(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, Field& field) {

		assert_true(pos.dominatedBy(field.size())) << "Position " << pos << " is outside field of size " << field.size();

		switch(universeProperties.useCase) {

		case UseCase::Dipole:
		{
			field[pos].E = { 0.0, 0.0, 0.0 };
			auto B = field[pos].B;
			auto location = getLocation(pos, universeProperties);

			double fac1 = -B0 * pow(Re, 3.0) / pow(allscale::api::user::data::sumOfSquares(location), 2.5);
			B.x = 3.0 * location.x * location.z * fac1;
			B.y = 3.0 * location.y * location.z * fac1;
			B.z = (2.0 * location.z * location.z - location.x * location.x - location.y * location.y) * fac1;

			field[pos].B = B;
		}
		break;

		case UseCase::ParticleWave:
		{
			// TODO: to provide
			break;
		}

		case UseCase::Test:
		{
			// TODO: to provide
			break;
		}

		default:
			assert_not_implemented() << "The specified use case is not supported yet!";
		}
	}

	/**
	* Initial version of the Field Solver: compute fields E and B for the Boris mover
	*
	* Fields are computed with respect to each particle position
	*/
	void computeFields(const Particle& p, Vector3<double> &E, Vector3<double> &B, const UseCase useCase) {
		switch(useCase) {
			case UseCase::Dipole:
			{
				E = { 0, 0, 0 };
				double fac1 = -B0 * pow(Re, 3.0) / pow(sumOfSquares(p.position), 2.5);
				B.x = 3.0 * p.position.x * p.position.z * fac1;
				B.y = 3.0 * p.position.y * p.position.z * fac1;
				B.z = (2.0 * p.position.z * p.position.z - p.position.x * p.position.x - p.position.y * p.position.y) * fac1;
				return;
			}
			case UseCase::ParticleWave:
			{
				E.x = sin(2.0 * M_PI * p.position.x) * cos(2.0 * M_PI * p.position.y);
				E.y = p.position.x * (1.0 - p.position.x) * p.position.y * (1.0 - p.position.y);
				E.z = p.position.x * p.position.x + p.position.z * p.position.z;
				B.x = 0.0;
				B.y = cos(2.0 * M_PI * p.position.z);
				B.z = sin(2.0 * M_PI * p.position.x);
				return;
			}
			default:
				assert_not_implemented() << "Unknown use case " << useCase << " for computeFields";
		}
	}

	/**
	* Explicit Field Solver: Fields are computed using leapfrog algorithm
	*/
	void solveFieldLeapfrog(const UniverseProperties& /*universeProperties*/, const utils::Coordinate<3>& pos, Field& field) {

		assert_true(pos.dominatedBy(field.size())) << "Position " << pos << " is outside universe of size " << field.size();

		// 1. Compute current density J as sum of particles density times particles velocity

		// 2. Compute electric field E using leapfrog with the time step delta t

		// 3. Compute magnetic field B using leapfrog with the time step delta t, but starts on delta t / 2
		//    Compute also magnetic field B on the center of each cell as average of all nodes
	}

} // end namespace ipic3d
