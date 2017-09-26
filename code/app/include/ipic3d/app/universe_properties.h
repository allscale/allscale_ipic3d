#pragma once

#include "allscale/api/user/data/grid.h"

#include "ipic3d/app/vector.h"
#include "ipic3d/app/parameters.h"

#include <map>

namespace ipic3d {

	using coordinate_type = allscale::api::user::data::GridPoint<3>;

	namespace {
		std::string getUseCaseName(const UseCase& useCase) {
			switch(useCase) {
			case UseCase::Dipole: return "Dipole";
			case UseCase::ParticleWave: return "ParticleWave";
			case UseCase::Test: return "Test";
			}
			// adhere to fallback policy in parameter parsing
			return "ParticleWave";
		}
	}

    /**
	* Printing support for the use case enum.
	*/
	std::ostream& operator<<(std::ostream& out, const UseCase& useCase) {
		return out << getUseCaseName(useCase);
	}

	// TODO: these defaults should not be externally visible
	namespace {

		// TODO: check and rename all these fields

		// Earth parameters
		static const double Re = 6378137.0; 		// meter (Earth radius)
		static const double B0 = 3.07e-5; 			// Tesla
		// Other parameters
		static const double e = 1.602176565e-19; 	// Elementary charge (Coulomb)
		static const double m = 1.672621777e-27; 	// Proton mass (kg)
		static const double c = 299792458.0; 		// speed of light (m/s)
		static const double K = 1e7 * e;    		// kinetic energy in eV converted to Joules

	}

	/**
	 * Holds all properties common to/uniform in a universe that are needed throughout the simulation
	 *
	*/
	struct UniverseProperties {

		// the use case
		UseCase useCase;
		// The size of this universe
		coordinate_type size;
		// The width of cells
		Vector3<double> cellWidth;
		// The timestep
		double dt;
		// Speed of light
		double speedOfLight;
		// planet radius
		double planetRadius;
		// object center
		Vector3<double> objectCenter;
		// physical point of origin corresponds to the front lower left corner
		// 		by defaul it is {0, 0, 0}
		Vector3<double> origin; 
		// magnetic field
		Vector3<double> magneticField;
		int FieldOutputCycle;
		std::string outputFileBaseName;

	    UniverseProperties(const UseCase& useCase = UseCase::Dipole, const coordinate_type& size = {1, 1, 1}, const Vector3<double>& cellWidth = {1.0, 1.0, 1.0},
			const double dt = 1.0, const double speedOfLight = 1.0, const double planetRadius = 0.0, const Vector3<double>& objectCenter = { 0.0, 0.0, 0.0 }, const Vector3<double>& origin = { 0.0, 0.0, 0.0 }, 
			const Vector3<double>& magneticField = {0.0, 0.0, 0.0}, const int FieldOutputCycle = 1)
	        : useCase(useCase), size(size), cellWidth(cellWidth), dt(dt), speedOfLight(speedOfLight), planetRadius(planetRadius), objectCenter(objectCenter), origin(origin), magneticField(magneticField), FieldOutputCycle(FieldOutputCycle) {
		    assert_true(size.x > 0 && size.y > 0 && size.z > 0) << "Expected positive non-zero universe size, but got " << size;
		    assert_true(cellWidth.x > 0 && cellWidth.y > 0 && cellWidth.z > 0) << "Expected positive non-zero cell widths, but got " << cellWidth;
		    assert_lt(0, dt) << "Expected positive non-zero time step, but got " << dt;
		    assert_lt(0, speedOfLight) << "Expected positive non-zero speed of light, but got " << speedOfLight;
		    assert_le(0, planetRadius) << "Expected positive or zero object radius, but got " << planetRadius;
	    }

		UniverseProperties(const Parameters& params)
			: useCase(params.useCase),
			size({ params.ncells.x, params.ncells.y, params.ncells.z }),
			cellWidth({ params.dspace.x, params.dspace.y, params.dspace.z }),
			dt( params.dt ),
			speedOfLight( params.c ),
			planetRadius( params.planetRadius ),
			objectCenter({ params.objectCenter.x, params.objectCenter.y, params.objectCenter.z }),
			magneticField({ params.B0.x, params.B0.y, params.B0.z }),
			FieldOutputCycle ( params.FieldOutputCycle ) 
		{
			origin.x = params.objectCenter.x - params.ncells.x * params.dspace.x / 2.0;
			origin.y = params.objectCenter.y - params.ncells.y * params.dspace.y / 2.0;
			origin.z = params.objectCenter.z - params.ncells.z * params.dspace.z / 2.0;
		}

	    friend std::ostream& operator<<(std::ostream& out, const UniverseProperties& props) {
			out << "Universe properties:" << std::endl;
			out << "\tUse Case: " << props.useCase << std::endl;
			out << "\tSize: " << props.size << std::endl;
			out << "\tCell width: " << props.cellWidth << std::endl;
			out << "\tTimestep: " << props.dt << std::endl;
			out << "\tSpeed of light: " << props.speedOfLight << std::endl;
			out << "\tPlanet radius: " << props.planetRadius << std::endl;
			out << "\tObject center: " << props.objectCenter << std::endl;
			out << "\tOrigin of the domain: " << props.origin << std::endl;
			out << "\tMagnetic field: " << props.magneticField << std::endl;
			out << "\tField's output cycle: " << props.FieldOutputCycle<< std::endl;
			return out;
		}

	};

	Vector3<double> getOriginOfCell(const coordinate_type& pos, const UniverseProperties& properties) {
		// do not check for strict domination, as pos could also refer to a field position
		assert_true(pos.dominatedBy(properties.size)) << "Position " << pos << " is outside universe of size " << properties.size;
		Vector3<double> tempPos{ (double)pos.x, (double)pos.y, (double)pos.z };
		return properties.origin + elementwiseProduct(tempPos, properties.cellWidth);
	}

	Vector3<double> getLocationForFields(const coordinate_type& pos, const UniverseProperties& properties) {
		// do not check for strict domination, as pos could also refer to a field position
		assert_true(pos.dominatedBy(properties.size + coordinate_type(1))) << "Position " << pos << " is outside universe of size " << properties.size;
		Vector3<double> tempPos{ (double)pos.x, (double)pos.y, (double)pos.z };
		return properties.origin + elementwiseProduct(tempPos, properties.cellWidth);
	}

	Vector3<double> getCenterOfCell(const coordinate_type& pos, const UniverseProperties& properties) {
		assert_true(pos.strictlyDominatedBy(properties.size)) << "Position " << pos << " is outside universe of size " << properties.size;
		return properties.origin + getOriginOfCell(pos, properties) + properties.cellWidth / 2;
	}

}
