#pragma once

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/data/vector.h"

#include "ipic3d/app/init_properties.h"
#include "ipic3d/app/universe_properties.h"
#include "ipic3d/app/utils/points.h"

namespace ipic3d {

	struct FieldNode {
		Vector3<double> E;				// electric field components defined on nodes
		Vector3<double> B;				// magnetic field components defined on nodes
		Vector3<double> Bext;			// external magnetic field on nodes
	};

	struct BcFieldCell {
		Vector3<double> Bc;				// magnetic field components defined on central points between nodes TODO: to clarify this
										// we assume that they are define at center of each cell
	};

	struct DensityNode {
		Vector3<double> J;				// current density
	};

	struct DensityCell {
		double rho;						// charge density
	};


	using Field = allscale::api::user::data::Grid<FieldNode,3>;	// a 3D grid of field nodes

	using BcField = allscale::api::user::data::Grid<BcFieldCell,3>;	// a 3D grid of magnetic field cells defined on centers

	using DensityCells = allscale::api::user::data::Grid<DensityCell,3>;	// a 3D grid of density cells

	using DensityNodes = allscale::api::user::data::Grid<DensityNode,3>;	// a 3D grid of density nodes

	// declaration
	void interpN2C(const utils::Coordinate<3>& pos, const Field& fields, BcField& bcfields);

	// definition
	Field initFields(const InitProperties& initProperties, const UniverseProperties& universeProperties) {

		using namespace allscale::api::user;

		utils::Size<3> start = 1;
		// determine the field size (grid size + 1 in each dimension)
		utils::Size<3> fieldSize = universeProperties.size + coordinate_type(3); // two for the two extra boundary cells and one as fields are defined on nodes of the cells
		utils::Size<3> workingFieldSize = universeProperties.size + coordinate_type(2);

		// the 3-D force fields
		Field fields(fieldSize);

		auto driftVel = initProperties.driftVelocity;
		assert_false(driftVel.empty()) << "Expected a drift velocity vector of at least length 1";
		auto ebc = crossProduct(driftVel[0], initProperties.magneticFieldAmplitude) * -1;

		pfor(start, workingFieldSize, [&](const utils::Coordinate<3>& cur) {

			// initialize electrical field
			fields[cur].E = ebc;

			// initialize magnetic field
			fields[cur].B = initProperties.magneticFieldAmplitude;

			// -- add earth model --

			switch(universeProperties.useCase) {

				case UseCase::Dipole: {

					// radius of the planet
					double a = universeProperties.objectRadius;

					// Dipole's Center
					auto objectCenter = universeProperties.objectCenter;
					// Node coordinates
					auto location = getLocationForFields(cur, universeProperties);

					auto diff = location - objectCenter;

					double r2 = allscale::api::user::data::sumOfSquares(diff);

					// Compute dipolar field B_ext

					if (r2 > a*a) {
						auto fac1 =  -universeProperties.magneticField.z * pow(a, 3) / pow(r2, 2.5);
						fields[cur].Bext.x = 3.0 * diff.x * diff.z * fac1;
						fields[cur].Bext.y = 3.0 * diff.y * diff.z * fac1;
						fields[cur].Bext.z = (2.0 * diff.z * diff.z - diff.x * diff.x - diff.y * diff.y) * fac1;
					} else { // no field inside the planet
						fields[cur].Bext = { 0.0, 0.0, 0.0 };
					}

					break;
				}

				case UseCase::ParticleWave: {

					fields[cur].Bext = { 0, 0, 0 };

					break;
				}

				default:
						assert_not_implemented()
							<< "The specified use case is not supported yet!";
			}

		});

		// return the produced field
		return std::move(fields);
	}

	BcField initBcFields(const UniverseProperties& universeProperties, const Field& field) {

		using namespace allscale::api::user;

		utils::Size<3> start = 1;
		utils::Size<3> fieldSize = universeProperties.size + coordinate_type(2); // two extra boundary cells
		utils::Size<3> workingFieldSize = universeProperties.size + coordinate_type(1);

		// the 3-D force fields
		BcField bcfield(fieldSize);

		pfor(start, workingFieldSize, [&](const utils::Coordinate<3>& cur) {

			// init magnetic field at centers
			interpN2C(cur, field, bcfield);
		});

		// return the produced field
		return std::move(bcfield);
	}


	/**
 	* calculate curl on nodes, given a vector field defined on central points
 	*/
	void computeCurlB(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, const BcField& bcfield, Vector3<double> &curl){
		// extract magnetic field values from centers of the cells
		Vector3<double> Bs[2][2][2];
		for(int i=0; i<2; i++) {
			for(int j=0; j<2; j++) {
				for(int k=0; k<2; k++) {
					utils::Coordinate<3> cur({pos[0]-i,pos[1]-j,pos[2]-k});
					Bs[i][j][k] = bcfield[cur].Bc;
				}
			}
		}

		// compute curl of E on central points
		double compZDY, compYDZ;
		double compXDZ, compZDX;
		double compYDX, compXDY;

		// curl - X
		compZDY = .25 * (Bs[0][0][0].z - Bs[0][1][0].z + Bs[0][0][1].z - Bs[0][1][1].z + Bs[1][0][0].z - Bs[1][1][0].z + Bs[1][0][1].z - Bs[1][1][1].z) / universeProperties.cellWidth.y;
		compYDZ = .25 * (Bs[0][0][0].y - Bs[0][0][1].y + Bs[1][0][0].y - Bs[1][0][1].y + Bs[0][1][0].y - Bs[0][1][1].y + Bs[1][1][0].y - Bs[1][1][1].y) / universeProperties.cellWidth.z;
		curl.x = compZDY - compYDZ;

		// curl - Y
		compXDZ = .25 * (Bs[0][0][0].x - Bs[0][0][1].x + Bs[1][0][0].x - Bs[1][0][1].x + Bs[0][1][0].x - Bs[0][1][1].x + Bs[1][1][0].x - Bs[1][1][1].x) / universeProperties.cellWidth.z;
		compZDX = .25 * (Bs[0][0][0].z - Bs[1][0][0].z + Bs[0][0][1].z - Bs[1][0][1].z + Bs[0][1][0].z - Bs[1][1][0].z + Bs[0][1][1].z - Bs[1][1][1].z) / universeProperties.cellWidth.x;
		curl.y = compXDZ - compZDX;

		// curl - Z
		compYDX = .25 * (Bs[0][0][0].y - Bs[1][0][0].y + Bs[0][0][1].y - Bs[1][0][1].y + Bs[0][1][0].y - Bs[1][1][0].y + Bs[0][1][1].y - Bs[1][1][1].y) / universeProperties.cellWidth.x;
		compXDY = .25 * (Bs[0][0][0].x - Bs[0][1][0].x + Bs[0][0][1].x - Bs[0][1][1].x + Bs[1][0][0].x - Bs[1][1][0].x + Bs[1][0][1].x - Bs[1][1][1].x) / universeProperties.cellWidth.y;
		curl.z = compYDX - compXDY;
	}

	/**
 	* calculate curl on central points, given a vector field defined on nodes
 	*/
	void computeCurlE(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, const Field& field, Vector3<double> &curl) {
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
 	* interpolate on nodes from central points for the magnetic field
 	*/
	void interpC2N(const utils::Coordinate<3>& pos, const BcField& bcfields, Field& fields) {
		// extract magnetic field values from centers of the cells
		Vector3<double> Bc(0);
		for(int i=0; i<2; i++) {
			for(int j=0; j<2; j++) {
				for(int k=0; k<2; k++) {
					utils::Coordinate<3> cur({pos[0]-i,pos[1]-j,pos[2]-k});
					Bc += bcfields[cur].Bc;
				}
			}
		}

		fields[pos].B = .125 * Bc;
	}

	/**
 	* interpolate on central points from nodes
 	*/
	void interpN2C(const utils::Coordinate<3>& pos, const Field& fields, BcField& bcfields) {
		// extract magnetic field values from nodes
		Vector3<double> Bn(0);
		for(int i=0; i<2; i++) {
			for(int j=0; j<2; j++) {
				for(int k=0; k<2; k++) {
					utils::Coordinate<3> cur({pos[0]+i,pos[1]+j,pos[2]+k});
					Bn += fields[cur].B;
				}
			}
		}

		bcfields[pos].Bc = .125 * Bn;
	}

	/**
	* Static filed solver in this case works as a push for simulation at the its beginning.
	* So that, fields are computed only once and then updated via interpolation
	* Computations should be done at centers
	*/
	void solveFieldStatically(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, Field& field) {

		assert_true(pos.dominatedBy(field.size())) << "Position " << pos << " is outside field of size " << field.size();

		switch(universeProperties.useCase) {

			case UseCase::Dipole:
			{
				field[pos].E = { 0.0, 0.0, 0.0 };
				auto B = field[pos].B;

				// Node coordinates
				auto location = getCenterOfCell(pos, universeProperties);

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
	* Explicit Field Solver: Fields are computed using forward approximation
	*/
	void solveFieldForward(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, DensityNodes& density, Field& field, BcField& bcfield) {

		assert_true(pos.dominatedBy(field.size())) << "Position " << pos << " is outside universe of size " << field.size();

		// 1. Compute current density J as sum of particles density times particles velocity

		// 2. Compute electric field E using leapfrog with the time step delta t

		// 3. Compute magnetic field B using leapfrog with the time step delta t, but starts on delta t / 2
		//    Compute also magnetic field B on the center of each cell as average of all nodes

		switch(universeProperties.useCase) {

			case UseCase::Dipole:
			{
				// 1. Compute E
				// 		curl of B
				Vector3<double> curlB;
				computeCurlB(universeProperties, pos, bcfield, curlB);

				// 		scale Jh by -4PI/c
				// 		sum curl B and Jh
				// 		scale the sum by dt
				// 		update E_{n+1} with the computed value
				field[pos].E += (curlB + density[pos - utils::Coordinate<3>(1)].J) * universeProperties.dt; // density needs to be shifted as pos corresponds to the fields position with a shift of one
 
				// 		TODO:Boundary conditions: periodic?
				//		periodic boundary conditions should be automatically supported as we added an extra row of cells around the grid


				// 2. Compute B
				// 		curl of E
				Vector3<double> curlE;
				computeCurlE(universeProperties, pos, field, curlE);

				//		scale curl by -c*dt
				//		TODO: check the speed of light
				//		update B_{n+1} on the center with the computed value
				bcfield[pos].Bc -= curlE * universeProperties.dt;

				//		TODO: Boundary conditions: periodic?
				//		periodic boundary conditions should be automatically supported as we added an extra row of cells around the grid

				// 		interpolate B from center to nodes
				interpC2N(pos, bcfield, field);

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

