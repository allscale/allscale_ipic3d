#pragma once

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/operator/pfor.h"

#include "ipic3d/app/vector.h"
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
		double rho;						// charge density
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

		// the 3D force fields
		Field fields(fieldSize);

		// the 3D density fields
		DensityNodes densityNodes(fieldSize);

		switch(universeProperties.useCase) {

			case UseCase::Dipole: {

				auto driftVel = initProperties.driftVelocity;
				assert_false(driftVel.empty()) << "Expected a drift velocity vector of at least length 1";
				auto ebc = -1.0 * crossProduct(driftVel[0], initProperties.magneticFieldAmplitude);

				// radius of the planet
				double a = universeProperties.planetRadius;

				// Dipole's Center
				auto objectCenter = universeProperties.objectCenter;

				pfor(start, workingFieldSize, [&](const utils::Coordinate<3>& cur) {

					// initialize rhos
					densityNodes[cur].rho = initProperties.rhoInit / (4.0 * M_PI); 

					// initialize electrical field
					fields[cur].E = ebc;

					// initialize magnetic field
					fields[cur].B = initProperties.magneticFieldAmplitude;

					// -- add earth model --

					// Node coordinates
					// TODO: double check cur - start
					// pos-start due to the fact that we have a ghost field
					auto location = getLocationForFields(cur-start, universeProperties);

					auto diff = location - objectCenter;

					double r2 = allscale::utils::sumOfSquares(diff);

					// Compute dipolar field B_ext
					if (r2 > a*a) {
						auto fac1 =  -universeProperties.magneticField.z * pow(a, 3) / pow(r2, 2.5);
						fields[cur].Bext.x = 3.0 * diff.x * diff.z * fac1;
						fields[cur].Bext.y = 3.0 * diff.y * diff.z * fac1;
						fields[cur].Bext.z = (2.0 * diff.z * diff.z - diff.x * diff.x - diff.y * diff.y) * fac1;
					} else { // no field inside the planet
						fields[cur].Bext = { 0.0, 0.0, 0.0 };
					}

					// TODO: investigate this (commented out in the original code)
					fields[cur].B += fields[cur].Bext;

				});

				break;
			}

			case UseCase::ParticleWave: {

				break;
			}

			default:
					assert_not_implemented()
						<< "The specified use case is not supported yet!";
		}

		// return the produced field
		return fields;
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
		return bcfield;
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
				// no need to update as Bext is computed only once during the initialization
				// the same holds for E and B

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
				field[pos].E += (universeProperties.speedOfLight * curlB - density[pos - utils::Coordinate<3>(1)].J) * universeProperties.dt; // density needs to be shifted as pos corresponds to the fields position with a shift of one
//				field[pos].E.x += (universeProperties.speedOfLight * curlB.z - density[pos - utils::Coordinate<3>(1)].J.x) * universeProperties.dt;
//				field[pos].E.y += (universeProperties.speedOfLight * curlB.z - density[pos - utils::Coordinate<3>(1)].J.y) * universeProperties.dt;
//				field[pos].E.y += (universeProperties.speedOfLight * (curlB.y - curlB.x) - density[pos - utils::Coordinate<3>(1)].J.z) * universeProperties.dt;

				// 2. Compute B
				// 		curl of E
				Vector3<double> curlE;
				computeCurlE(universeProperties, pos, field, curlE);

				//		scale curl by -c*dt
				//		update B_{n+1} on the center with the computed value
				bcfield[pos].Bc -= universeProperties.speedOfLight * curlE * universeProperties.dt;
//				field[pos].B.x -= universeProperties.speedOfLight * curlE.z * universeProperties.dt;
//				field[pos].B.y -= universeProperties.speedOfLight * curlE.z * universeProperties.dt;
//				field[pos].B.z -= universeProperties.speedOfLight * (curlE.y - curlB.x) * universeProperties.dt;

				// 		Boundary conditions: periodic are supported automatically supported as we added an extra row of cells around the grid

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
	void solveFieldLeapfrog(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, DensityNodes& density, Field& field, BcField& bcfield) {

		assert_true(pos.dominatedBy(field.size())) << "Position " << pos << " is outside universe of size " << field.size();

		// compute electric field E using leapfrog with the time step delta t

		// compute magnetic field B using leapfrog with the time step delta t, but starts on delta t / 2
		// compute also magnetic field B on the center of each cell as average of all nodes

		switch(universeProperties.useCase) {

			case UseCase::Dipole:
			{
				//	Compute transverse magnetic (TM) sets
				bcfield[pos].Bc.z = bcfield[pos].Bc.z - universeProperties.speedOfLight * universeProperties.dt *( (field[pos+utils::Coordinate<3>({1,0,0})].E.y - field[pos].E.y) / universeProperties.cellWidth.x + (field[pos+utils::Coordinate<3>({0,1,0})].E.x - field[pos].E.x) / universeProperties.cellWidth.y );

				field[pos].E.x = field[pos].E.x + universeProperties.speedOfLight * universeProperties.dt * (bcfield[pos].Bc.z - bcfield[pos+utils::Coordinate<3>({0,-1,0})].Bc.z) / universeProperties.cellWidth.x - universeProperties.dt * density[pos].J.x; 

				field[pos].E.y = field[pos].E.y - universeProperties.speedOfLight * universeProperties.dt * (bcfield[pos].Bc.z - bcfield[pos+utils::Coordinate<3>({-1,0,0})].Bc.z) / universeProperties.cellWidth.y - universeProperties.dt * density[pos].J.y; 


				//	Compute transverse electric (TE) sets
				field[pos].E.z = field[pos].E.z + universeProperties.dt * ( universeProperties.speedOfLight * (bcfield[pos+utils::Coordinate<3>({1,0,0})].Bc.y- bcfield[pos].Bc.y) / universeProperties.cellWidth.x - universeProperties.speedOfLight * (bcfield[pos+utils::Coordinate<3>({0,1,0})].Bc.x- bcfield[pos].Bc.x) / universeProperties.cellWidth.y - density[pos].J.z );

				bcfield[pos].Bc.x = bcfield[pos].Bc.x - universeProperties.speedOfLight * universeProperties.dt * (field[pos].E.z - field[pos+utils::Coordinate<3>({0,-1,0})].E.z) / universeProperties.cellWidth.y; 

				bcfield[pos].Bc.y = bcfield[pos].Bc.y - universeProperties.speedOfLight * universeProperties.dt * (field[pos].E.z - field[pos+utils::Coordinate<3>({-1,0,0})].E.z) / universeProperties.cellWidth.x; 

				//	Boundary conditions: periodic are supported automatically supported as we added an extra row of cells around the grid

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
 	* Populate and update fields values on boundaries
 	*/
	void updateFieldsOnBoundaries(Field& field, BcField& bcfield) {
		// update on the x = 0 and x = N face
		int fieldSize = (int)field.size().y;
		int bcfieldSize = (int)bcfield.size().y;
		for (int i = 1; i < fieldSize - 1; i++) { // for y
			for (int j = 1; j < fieldSize - 1; j++) { // for z
				utils::Coordinate<3> pos0{0, i, j};
				utils::Coordinate<3> pos1{1, i, j};
				utils::Coordinate<3> posN1{fieldSize-2, i, j};
				utils::Coordinate<3> posN{fieldSize-1, i, j};
				
				field[pos0] = field[posN1];
				field[posN] = field[pos1];
			}
		}
		for (int i = 1; i < bcfieldSize - 1; i++) { // for y
			for (int j = 1; j < bcfieldSize - 1; j++) { // for z
				utils::Coordinate<3> pos0{0, i, j};
				utils::Coordinate<3> pos1{1, i, j};
				utils::Coordinate<3> posN1{bcfieldSize-2, i, j};
				utils::Coordinate<3> posN{bcfieldSize-1, i, j};

				bcfield[pos0] = bcfield[posN1];
				bcfield[posN] = bcfield[pos1];
			}
		}

		// update on the y = 0 face
		for (int i = 1; i < fieldSize - 1; i++) { // for x
			for (int j = 1; j < fieldSize - 1; j++) { // for z
				utils::Coordinate<3> pos0{i, 0, j};
				utils::Coordinate<3> pos1{i, 1, j};
				utils::Coordinate<3> posN1{i, fieldSize-2, j};
				utils::Coordinate<3> posN{i, fieldSize-1, j};
				
				field[pos0] = field[posN1];
				field[posN] = field[pos1];
			}
		}
		for (int i = 1; i < bcfieldSize - 1; i++) { // for y
			for (int j = 1; j < bcfieldSize - 1; j++) { // for z
				utils::Coordinate<3> pos0{i, 0, j};
				utils::Coordinate<3> pos1{i, 1, j};
				utils::Coordinate<3> posN1{i, bcfieldSize-2, j};
				utils::Coordinate<3> posN{i, bcfieldSize-1, j};

				bcfield[pos0] = bcfield[posN1];
				bcfield[posN] = bcfield[pos1];
			}
		}

		// update on the z = 0 face
		for (int i = 1; i < fieldSize - 1; i++) { // for y
			for (int j = 1; j < fieldSize - 1; j++) { // for z
				utils::Coordinate<3> pos0{i, j, 0};
				utils::Coordinate<3> pos1{i, j, 1};
				utils::Coordinate<3> posN1{i, j, fieldSize-2};
				utils::Coordinate<3> posN{i, j, fieldSize-1};
				
				field[pos0] = field[posN1];
				field[posN] = field[pos1];
			}
		}
		for (int i = 1; i < bcfieldSize - 1; i++) { // for y
			for (int j = 1; j < bcfieldSize - 1; j++) { // for z
				utils::Coordinate<3> pos0{i, j, 0};
				utils::Coordinate<3> pos1{i, j, 1};
				utils::Coordinate<3> posN1{i, j, bcfieldSize-2};
				utils::Coordinate<3> posN{i, j, bcfieldSize-1};

				bcfield[pos0] = bcfield[posN1];
				bcfield[posN] = bcfield[pos1];
			}
		}
	}

	// compute the electric field energy
	double getEenergy(const Field& field, const UniverseProperties& universeProperties){
		auto fieldSize = field.size();
		auto fieldStart = utils::Coordinate<3>(1);
		auto fieldEnd = fieldSize - utils::Coordinate<3>(1); // one because a shift due to the boundary conditions
		
		double sum = 0.0;
		double vol = 0.5 * universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;
		double fourPI = 4.0 * M_PI;
		allscale::api::user::pfor(fieldStart, fieldEnd, [&](const utils::Coordinate<3>& pos){
			auto e = field[pos].E;
			sum += vol * ( e.x * e.x + e.y * e.y + e.z * e.z ) / (fourPI);
		});

		return sum;
	} 

	// compute the magnetic field energy
	double getBenergy(const Field& field, const UniverseProperties& universeProperties){
		auto fieldSize = field.size();
		auto fieldStart = utils::Coordinate<3>(1);
		auto fieldEnd = fieldSize - utils::Coordinate<3>(1); // one because a shift due to the boundary conditions
		
		double sum = 0.0;
		double vol = 0.5 * universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;
		double fourPI = 4.0 * M_PI;
		allscale::api::user::pfor(fieldStart, fieldEnd, [&](const utils::Coordinate<3>& pos){
			auto b = field[pos].B;
			sum += vol * ( b.x * b.x + b.y * b.y + b.z * b.z ) / (fourPI);
		});

		return sum;
	} 

} // end namespace ipic3d

