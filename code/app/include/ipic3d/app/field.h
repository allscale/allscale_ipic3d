#pragma once

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/data/vector.h"

#include "ipic3d/app/universe_properties.h"
#include "ipic3d/app/utils/points.h"

namespace ipic3d {

	struct FieldNode {
		Vector3<double> E;				// electric field components defined on nodes
		Vector3<double> B;				// magnetic field components defined on nodes
		Vector3<double> Bc;				// magnetic field components defined on central points between nodes TODO: to clarify this
		Vector3<double> Bext;			// external magnetic field on nodes
	};


	using Field = allscale::api::user::data::Grid<FieldNode,3>;	// a 3D grid of field nodes

	/**
 	* calculate curl on nodes, given a vector field defined on central points
 	*/
	void curlC2N(){

	}

	/**
 	* calculate curl on central points, given a vector field defined on nodes
 	*/
	void curlN2C(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, const Field& field, Vector3<double> &curl) {
		// extract electric field values from nodes of the cell
		Vector3<double> Es[2][2][2];
		for(int i=0; i<2; i++) {
			for(int j=0; j<2; j++) {
				for(int k=0; k<2; k++) {
					utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});
					Es[i][j][k] = field[cur].E;
				}
			}
		}

		// compute curl of E on central points
		curl = 0;
		double compZDY, compYDZ;
		double compXDZ, compZDX;
		double compYDX, compXDY;

		// curl - X
		compZDY = .25 * (Es[0][1][0].z - Es[0][0][0].z + Es[0][1][1].z - Es[0][0][1].z + Es[1][1][0].z - Es[1][0][0].z + Es[1][1][1].z - Es[1][0][1].z) / universeProperties.cellWidth.y;
		compYDZ = .25 * (Es[0][0][1].y - Es[0][0][0].y + Es[1][0][1].y - Es[1][0][0].y + Es[0][1][1].y - Es[0][1][0].y + Es[1][1][1].y - Es[1][1][0].y) / universeProperties.cellWidth.z;
		curl.x = compZDY - compYDZ; 

		// curl - Y
		compXDZ = .25 * (Es[0][0][1].x - Es[0][0][0].x + Es[1][0][1].x - Es[1][0][0].x + Es[0][1][1].x - Es[0][1][0].x + Es[1][1][1].x - Es[1][1][0].x) / universeProperties.cellWidth.z;
		compZDX = .25 * (Es[1][0][0].z - Es[0][0][0].z + Es[1][0][1].z - Es[0][0][1].z + Es[1][1][0].z - Es[0][1][0].z + Es[1][1][1].z - Es[0][1][1].z) / universeProperties.cellWidth.x;
		curl.y = compXDZ - compZDX; 

		// curl - Z
		compYDX = .25 * (Es[1][0][0].y - Es[0][0][0].y + Es[1][0][1].y - Es[0][0][1].y + Es[1][1][0].y - Es[0][1][0].y + Es[1][1][1].y - Es[0][1][1].y) / universeProperties.cellWidth.x;
		compXDY = .25 * (Es[0][1][0].x - Es[0][0][0].x + Es[0][1][1].x - Es[0][0][1].x + Es[1][1][0].x - Es[1][0][0].x + Es[1][1][1].x - Es[1][0][1].x) / universeProperties.cellWidth.y;
		curl.z = compYDX - compXDY; 
	}

	/**
	* TODO: provide a description of how the static field solver works
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
				break;
			}

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
	* Explicit Field Solver: Fields are computed using forward approximation
	*/
	void solveFieldForward(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, Field& field) {

		assert_true(pos.dominatedBy(field.size())) << "Position " << pos << " is outside universe of size " << field.size();

		// 1. Compute current density J as sum of particles density times particles velocity

		// 2. Compute electric field E using leapfrog with the time step delta t

		// 3. Compute magnetic field B using leapfrog with the time step delta t, but starts on delta t / 2
		//    Compute also magnetic field B on the center of each cell as average of all nodes

		switch(universeProperties.useCase) {

			case UseCase::Dipole:
			{
				// 1. Compute E
				// 		curlC2N()
				// 		scale Jh by -4PI/c
				// 		sum curl B and Jh
				// 		scale the sum by dt
				// 		update E_{n+1} with the computed value
				// 		Boundary conditions: periodic?
				
				// 2. Compute B
				// 		curlN2C()
				Vector3<double> curlE;
				curlN2C(universeProperties, pos, field, curlE);

				//		scale curl by -c*dt
				//		TODO: check the speed of light
				curlE -= curlE * universeProperties.dt;

				//		update B_{n+1} on the center with the computed value
				//		TODO: B should be defined in the center
				field[pos].Bc = field[pos].Bc + curlE;

				//		TODO: Boundary conditions: periodic?
						
				//		interpC2N
				// 		interpolate B from center to nodes			
				break;
			}

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
	* Explicit Field Solver: Fields are computed using leapfrog algorithm
	*/
	void solveFieldLeapfrog(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, Field& field) {

		assert_true(pos.dominatedBy(field.size())) << "Position " << pos << " is outside universe of size " << field.size();

		// 1. Compute current density J as sum of particles density times particles velocity

		// 2. Compute electric field E using leapfrog with the time step delta t

		// 3. Compute magnetic field B using leapfrog with the time step delta t, but starts on delta t / 2
		//    Compute also magnetic field B on the center of each cell as average of all nodes

		switch(universeProperties.useCase) {

			case UseCase::Dipole:
			{
				// TODO
				break;
			}

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

} // end namespace ipic3d

