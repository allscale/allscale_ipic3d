#pragma once

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/user/algorithm/preduce.h"

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
	};

	//struct DensityCell {
	//	Vector3<double> rho;			// charge density
	//};


	using Field = allscale::api::user::data::Grid<FieldNode,3>;	// a 3D grid of field nodes

	using BcField = allscale::api::user::data::Grid<BcFieldCell,3>;	// a 3D grid of magnetic field cells defined on centers

	using CurrentDensity = allscale::api::user::data::Grid<DensityNode,3>;	// a 3D grid of density nodes

	// declaration
	void interpN2C(const utils::Coordinate<3>& pos, const Field& fields, BcField& bcfields);

	FieldNode getDipoleFieldAt(const Vector3<double>& location, const InitProperties& initProperties, const UniverseProperties& universeProperties) {

		auto driftVel = initProperties.driftVelocity;
		assert_false(driftVel.empty()) << "Expected a drift velocity vector of at least length 1";
		auto ebc = -1.0 * crossProduct(driftVel[0], initProperties.magneticField);

		FieldNode res;

		// initialize electrical field
		res.E = ebc;

		// -- add earth model --

		// radius of the planet
		double a = universeProperties.planetRadius;

		auto diff = location - universeProperties.objectCenter;

		double r2 = allscale::utils::sumOfSquares(diff);

		// Compute dipolar field B_ext
		if (r2 > a*a) {
			auto fac1 =  -universeProperties.externalMagneticField.z * pow(a, 3) / pow(r2, 2.5);
			res.Bext.x = 3.0 * diff.x * diff.z * fac1;
			res.Bext.y = 3.0 * diff.y * diff.z * fac1;
			res.Bext.z = (2.0 * diff.z * diff.z - diff.x * diff.x - diff.y * diff.y) * fac1;
		} else { // no field inside the planet
			res.Bext = { 0.0, 0.0, 0.0 };
		}

		// initialize magnetic field
		res.B = initProperties.magneticField + res.Bext;

		// done
		return res;
	}

	// definition
	Field initFields(const InitProperties& initProperties, const UniverseProperties& universeProperties) {

		using namespace allscale::api::user::algorithm;

		utils::Size<3> start = 1;
		// determine the field size (grid size + 1 in each dimension)
		utils::Size<3> fieldSize = universeProperties.size + coordinate_type(3); // two for the two extra boundary cells and one as fields are defined on nodes of the cells
		utils::Size<3> workingFieldSize = universeProperties.size + coordinate_type(2);

		// the 3D force fields
		Field fields(fieldSize);

		switch(universeProperties.useCase) {

			case UseCase::Dipole: {

				pfor(start, workingFieldSize, [=,&fields](const utils::Coordinate<3>& cur) {

					// Node coordinates
					// pos-start due to the fact that we have a ghost field
					auto location = getLocationForFields(cur-start, universeProperties);

					// init current field cell
					fields[cur] = getDipoleFieldAt(location, initProperties, universeProperties);

				});

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

		using namespace allscale::api::user::algorithm;

		utils::Size<3> start = 1;
		utils::Size<3> fieldSize = universeProperties.size + coordinate_type(2); // two extra boundary cells
		utils::Size<3> workingFieldSize = universeProperties.size + coordinate_type(1);

		// the 3-D force fields
		BcField bcfield(fieldSize);

		pfor(start, workingFieldSize, [=,&field,&bcfield](const utils::Coordinate<3>& cur) {
			// init magnetic field at centers
			interpN2C(cur, field, bcfield);
		});

		// return the produced field
		return bcfield;
	}

	CurrentDensity initCurrentDensity(const UniverseProperties& universeProperties) {

		using namespace allscale::api::user::algorithm;

		utils::Size<3> start = 0;
		utils::Size<3> densitySize = universeProperties.size + coordinate_type(1);

		// the 3D current density
		CurrentDensity currentDensity(densitySize);

		pfor(start, densitySize, [=,&currentDensity](const utils::Coordinate<3>& cur) {

			// initialize current density on nodes
			currentDensity[cur].J = { 0.0, 0.0, 0.0 };
		});

		return currentDensity;
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

			default:
				assert_not_implemented() << "The specified use case is not supported yet!";
		}
	}

	/**
	* Explicit Field Solver: Fields are computed using forward approximation
	*/
	void solveFieldForward(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, CurrentDensity& density, Field& field, BcField& bcfield) {

		assert_true(pos.dominatedBy(field.size())) << "Position " << pos << " is outside universe of size " << field.size();

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

				// do this check to not-udpate the ghost cells here
				// TODO: a better idea to avoid this check
				if (pos < (bcfield.size() - utils::Coordinate<3>(1))) {
					// 2. Compute B
					// 		curl of E
					Vector3<double> curlE;
					computeCurlE(universeProperties, pos, field, curlE);

					//		scale curl by -c*dt
					//		update B_{n+1} on the center with the computed value
					bcfield[pos].Bc -= universeProperties.speedOfLight * curlE * universeProperties.dt;
				}

				// 		Boundary conditions: periodic are supported automatically supported as we added an extra row of cells around the grid

				// 		interpolate B from center to nodes
				interpC2N(pos, bcfield, field);

				break;
			}

			default:
				assert_not_implemented() << "The specified use case is not supported yet!";
		}
	}

	/**
	* Explicit Field Solver: Fields are computed using leapfrog algorithm
	*/
	void solveFieldLeapfrog(const UniverseProperties& universeProperties, const utils::Coordinate<3>& pos, CurrentDensity& density, Field& field, BcField& bcfield) {

		assert_true(pos.dominatedBy(field.size())) << "Position " << pos << " is outside universe of size " << field.size();

		// 1. Compute electric field E using leapfrog with the time step delta t

		// 2. Compute magnetic field B using leapfrog with the time step delta t, but starts on delta t / 2
		//    Compute also magnetic field B on the center of each cell as average of all nodes

		switch(universeProperties.useCase) {

			case UseCase::Dipole:
			// new approach to appear here
			{
				if (1) {
					// this approach is adopted from the article 'The Plasma Simulation Code: A modern particle-in-cell code with load-balancing and GPU support' by K. Germaschewski et al.
					// Compute E field
					field[pos].E.x = field[pos].E.x + universeProperties.dt * ( universeProperties.speedOfLight * universeProperties.speedOfLight * ( (bcfield[pos].Bc.z - bcfield[pos-utils::Coordinate<3>({0,1,0})].Bc.z) / universeProperties.cellWidth.y - (bcfield[pos].Bc.y - bcfield[pos-utils::Coordinate<3>({0,0,1})].Bc.y) / universeProperties.cellWidth.z) - density[pos - utils::Coordinate<3>(1)].J.x ); 

					field[pos].E.y = field[pos].E.y + universeProperties.dt * ( universeProperties.speedOfLight * universeProperties.speedOfLight * ( (bcfield[pos].Bc.x - bcfield[pos-utils::Coordinate<3>({0,0,1})].Bc.x) / universeProperties.cellWidth.z - (bcfield[pos].Bc.z - bcfield[pos-utils::Coordinate<3>({1,0,0})].Bc.z) / universeProperties.cellWidth.x) - density[pos - utils::Coordinate<3>(1)].J.y ); 

					field[pos].E.z = field[pos].E.z + universeProperties.dt * ( universeProperties.speedOfLight * universeProperties.speedOfLight * ( (bcfield[pos].Bc.y - bcfield[pos-utils::Coordinate<3>({1,0,0})].Bc.y) / universeProperties.cellWidth.x - (bcfield[pos].Bc.x - bcfield[pos-utils::Coordinate<3>({0,1,0})].Bc.x) / universeProperties.cellWidth.y) - density[pos - utils::Coordinate<3>(1)].J.z ); 

					//	Compute B field
					if (pos < (bcfield.size() - utils::Coordinate<3>(1))) {
						bcfield[pos].Bc.x = bcfield[pos].Bc.x - universeProperties.dt * ( (field[pos+utils::Coordinate<3>({0,1,0})].E.z - field[pos].E.z) / universeProperties.cellWidth.y - (field[pos+utils::Coordinate<3>({0,0,1})].E.y - field[pos].E.y) / universeProperties.cellWidth.z );

						bcfield[pos].Bc.y = bcfield[pos].Bc.y - universeProperties.dt * ( (field[pos+utils::Coordinate<3>({0,0,1})].E.x - field[pos].E.x) / universeProperties.cellWidth.z - (field[pos+utils::Coordinate<3>({1,0,0})].E.z - field[pos].E.z) / universeProperties.cellWidth.x );

						bcfield[pos].Bc.z = bcfield[pos].Bc.z - universeProperties.dt * ( (field[pos+utils::Coordinate<3>({1,0,0})].E.y - field[pos].E.y) / universeProperties.cellWidth.x - (field[pos+utils::Coordinate<3>({0,1,0})].E.x - field[pos].E.x) / universeProperties.cellWidth.y );
					}

				} else {
					// this approach is from the Birdsall book 'Plasma Physics via Computer Simulation'
					//	Compute transverse magnetic (TM) sets
					if (pos < (bcfield.size() - utils::Coordinate<3>(1))) {
						bcfield[pos].Bc.z = bcfield[pos].Bc.z - universeProperties.speedOfLight * universeProperties.dt *( (field[pos+utils::Coordinate<3>({1,0,0})].E.y - field[pos].E.y) / universeProperties.cellWidth.x + (field[pos+utils::Coordinate<3>({0,1,0})].E.x - field[pos].E.x) / universeProperties.cellWidth.y );
					}

					field[pos].E.x = field[pos].E.x + universeProperties.speedOfLight * universeProperties.dt * (bcfield[pos].Bc.z - bcfield[pos+utils::Coordinate<3>({0,-1,0})].Bc.z) / universeProperties.cellWidth.x - universeProperties.dt * density[pos - utils::Coordinate<3>(1)].J.x; 

					field[pos].E.y = field[pos].E.y - universeProperties.speedOfLight * universeProperties.dt * (bcfield[pos].Bc.z - bcfield[pos+utils::Coordinate<3>({-1,0,0})].Bc.z) / universeProperties.cellWidth.y - universeProperties.dt * density[pos - utils::Coordinate<3>(1)].J.y; 


					//	Compute transverse electric (TE) sets
					field[pos].E.z = field[pos].E.z + universeProperties.dt * ( universeProperties.speedOfLight * (bcfield[pos+utils::Coordinate<3>({1,0,0})].Bc.y- bcfield[pos].Bc.y) / universeProperties.cellWidth.x - universeProperties.speedOfLight * (bcfield[pos+utils::Coordinate<3>({0,1,0})].Bc.x- bcfield[pos].Bc.x) / universeProperties.cellWidth.y - density[pos - utils::Coordinate<3>(1)].J.z );

					if (pos < (bcfield.size() - utils::Coordinate<3>(1))) {
						bcfield[pos].Bc.x = bcfield[pos].Bc.x - universeProperties.speedOfLight * universeProperties.dt * (field[pos].E.z - field[pos+utils::Coordinate<3>({0,-1,0})].E.z) / universeProperties.cellWidth.y; 

						bcfield[pos].Bc.y = bcfield[pos].Bc.y - universeProperties.speedOfLight * universeProperties.dt * (field[pos].E.z - field[pos+utils::Coordinate<3>({-1,0,0})].E.z) / universeProperties.cellWidth.x; 
					}
				}

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

		auto update = [](const auto& index, const auto& end, auto& field) {

			// update on the x = 0 face
			field[{ 0, index.x, index.y }] = field[{ end - 1, index.x, index.y }];
			field[{ end, index.x, index.y }] = field[{ 1, index.x, index.y }];

			// update on the y = 0 face
			field[{ index.x, 0, index.y }] = field[{ index.x, end - 1, index.y }];
			field[{ index.x, end, index.y }] = field[{ index.x, 1, index.y }];

			// update on the z = 0 face
			field[{ index.x, index.y, 0 }] = field[{ index.x, index.y, end - 1 }];
			field[{ index.x, index.y, end }] = field[{ index.x, index.y, 1 }];
		};

		// Note: Assumes that the field sizes are equal in each dimension
		int fieldEnd = (int)field.size().x - 1;
		int bcfieldEnd = (int)bcfield.size().x - 1;
		allscale::utils::Vector<int, 2> fullField(fieldEnd);
		allscale::utils::Vector<int, 2> fullBcField(bcfieldEnd);
		allscale::utils::Vector<int, 2> start(1);

		// TODO: parallel version pending
		allscale::api::user::algorithm::detail::forEach(start, fullField, [&](const auto& index) {
			update(index, fieldEnd, field);
		});

		allscale::api::user::algorithm::detail::forEach(start, fullBcField, [&](const auto& index) {
			update(index, bcfieldEnd, bcfield);
		});

	}

	// compute the electric field energy
	template<typename Accessor>
	double getFieldEnergy(const Field& field, const UniverseProperties& universeProperties, const Accessor& accessor){
		auto fieldStart = utils::Coordinate<3>(1);
		auto fieldEnd = field.size() - utils::Coordinate<3>(1); // one because a shift due to the boundary conditions
		
		const double vol = 0.5 * universeProperties.cellWidth.x * universeProperties.cellWidth.y * universeProperties.cellWidth.z;
		const double fourPI = 16.0 * atan(1.0);

		auto map = [&](const coordinate_type& index, double& res) {
			res += vol * allscale::utils::sumOfSquares(accessor(field, index)) / (fourPI);
		};

		auto reduce = [&](const double& a, const double& b) { return a + b; };
		auto init = []() { return 0.0; };

		return allscale::api::user::algorithm::preduce(fieldStart, fieldEnd, map, reduce, init).get();
	} 

	/**
	* This function outputs all field values
	*/
	void outputFieldGrids(const Field& field, const BcField& bcField, const std::string& outputFilename) {
		// TODO: implement output facilities for large problems
		assert_le(field.size(), (coordinate_type{ 32,32,32 })) << "Unable to dump data for such a large field at this time";

		// output field values
		allscale::api::user::algorithm::async([=, &field]() {
			auto& manager = allscale::api::core::FileIOManager::getInstance();
			auto text = manager.createEntry(outputFilename);
			auto out = manager.openOutputStream(text);

			// output dimensions
			out << field.size() << "\n";

			for(std::int64_t i = 0; i < field.size().x; ++i) {
				for(std::int64_t j = 0; j < field.size().y; ++j) {
					for(std::int64_t k = 0; k < field.size().z; ++k) {
						coordinate_type p{ i,j,k };
						// write index
						out << p.x << "," << p.y << "," << p.z << ":";
						// write data
						out << field[p].E << "|" << field[p].B << "|" << field[p].Bext << "\n";
					}
				}
			}
			out << "\n";
			manager.close(out);
		}).wait();
		//allscale::api::user::algorithm::pfor(field.size(), [&](const auto& index) {
		//	streamObject.atomic([&](auto& out) {
		//		// write index
		//		out << index.x << "," << index.y << "," << index.z << ":";
		//		// write data
		//		out << field[index].E << "|" << field[index].B << "|" << field[index].Bext << "\n";
		//	});
		//});

		// output bc field values
		allscale::api::user::algorithm::async([=, &bcField]() {
			auto& manager = allscale::api::core::FileIOManager::getInstance();
			auto text = manager.createEntry(outputFilename);
			auto out = manager.openOutputStream(text);

			// output dimensions
			out << bcField.size() << "\n";

			for(std::int64_t i = 0; i < bcField.size().x; ++i) {
				for(std::int64_t j = 0; j < bcField.size().y; ++j) {
					for(std::int64_t k = 0; k < bcField.size().z; ++k) {
						coordinate_type p{ i,j,k };
						// write index
						out << p.x << "," << p.y << "," << p.z << ":";
						// write data
						out << bcField[p].Bc << "\n";
					}
				}
			}
			out << "\n";
			manager.close(out);
		}).wait();
		//allscale::api::user::algorithm::pfor(bcField.size(), [&](const auto& index) {
		//	streamObject.atomic([&](auto& out) {
		//		// write index
		//		out << index.x << "," << index.y << "," << index.z << ":";
		//		out << bcField[index].Bc << "\n"; 
		//	});
		//});
	}

} // end namespace ipic3d

