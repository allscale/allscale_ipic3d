#pragma once

#include "allscale/utils/printer/vectors.h"

#include "ipic3d/app/vector.h"
#include "ipic3d/app/parameters.h"

namespace ipic3d {

	/**
	 * Holds all properties that are needed for initialization only
	 *
	*/
	struct InitProperties {

		// the number of time steps
		unsigned numSteps;

		// the number of particles per cell per dimension, one entry per species
		std::vector<Vector3<unsigned>> particlesPerCell;

		// drift velocity
		std::vector<Vector3<double>> driftVelocity;

		// initial magnetic field
		Vector3<double> magneticField;

		// initial exteranal magnetic field
		Vector3<double> externalMagneticField;

		// charge density defined on nodes
		double rhoInit;

		InitProperties(const unsigned numSteps = 1, const std::vector<Vector3<unsigned>>& particlesPerCell = {},
			const std::vector<Vector3<double>>& driftVelocity = {}, const Vector3<double>& magneticField = { 0,0,0 }, const Vector3<double>& externalMagneticField = { 0,0,0 }, const double rhoInit = 1.0)
			: numSteps(numSteps), particlesPerCell(particlesPerCell), driftVelocity(driftVelocity), magneticField(magneticField), externalMagneticField(externalMagneticField), rhoInit(rhoInit) {}

		InitProperties(const Parameters& params) {

			numSteps = params.ncycles;

			for(int i = 0; i < params.ns; i++) {
				driftVelocity.push_back({ params.u0[i], params.v0[i], params.w0[i] });
			}

			for(int i = 0; i < (params.ns); i++) {
				particlesPerCell.push_back({ (unsigned)params.npcelx[i], (unsigned)params.npcely[i], (unsigned)params.npcelz[i] });
			}

			magneticField = { params.B0.x, params.B0.y, params.B0.z };
			externalMagneticField = { params.B1.x, params.B1.y, params.B1.z };

			rhoInit = params.rhoInit[0];
		}

	    friend std::ostream& operator<<(std::ostream& out, const InitProperties& props) {
			out << "InitProperties:" << std::endl;
			out << "\tNumber of time steps: " << props.numSteps << std::endl;
			out << "\tNumber of particles per cell: " << props.particlesPerCell << std::endl;
			out << "\tDrift velocity: " << props.driftVelocity << std::endl;
			out << "\tMagnetic fielde: " << props.magneticField << std::endl;
			out << "\tExternal magnetic fielde: " << props.externalMagneticField << std::endl;
			out << "\tCharge density: " << props.rhoInit << std::endl;
			return out;
		}

	};

}
