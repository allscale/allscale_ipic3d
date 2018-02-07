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
		std::uint64_t numSteps;

		// the number of particles per cell per dimension, one entry per species
		std::vector<Vector3<unsigned>> particlesPerCell;

		// drift velocity
		std::vector<Vector3<double>> driftVelocity;

		// magnetic field amplitude
		Vector3<double> magneticFieldAmplitude;

		// charge density defined on nodes
		double rhoInit;

		InitProperties(const std::uint64_t numSteps = 1, const std::vector<Vector3<unsigned>>& particlesPerCell = {},
			const std::vector<Vector3<double>>& driftVelocity = {}, const Vector3<double>& magneticFieldAmplitude = { 0,0,0 }, const double rhoInit = 1.0)
			: numSteps(numSteps), particlesPerCell(particlesPerCell), driftVelocity(driftVelocity), magneticFieldAmplitude(magneticFieldAmplitude), rhoInit(rhoInit) {}

		InitProperties(const Parameters& params) {

			numSteps = params.ncycles;

			for(int i = 0; i < params.ns; i++) {
				driftVelocity.push_back({ params.u0[i], params.v0[i], params.w0[i] });
			}

			for(int i = 0; i < (params.ns); i++) {
				assert_true(params.npcelx[i] >= 0 && params.npcely[i] >= 0 && params.npcelz[i] >= 0) << "Expected positive number of particles.";
				particlesPerCell.push_back({ (unsigned)params.npcelx[i], (unsigned)params.npcely[i], (unsigned)params.npcelz[i] });
			}

			magneticFieldAmplitude = { params.B0.x, params.B0.y, params.B0.z };

			rhoInit = params.rhoInit[0];
		}

	    friend std::ostream& operator<<(std::ostream& out, const InitProperties& props) {
			out << "InitProperties:" << std::endl;
			out << "\tNumber of time steps: " << props.numSteps << std::endl;
			out << "\tNumber of particles per cell: " << props.particlesPerCell << std::endl;
			out << "\tDrift velocity: " << props.driftVelocity << std::endl;
			out << "\tMagnetic field amplitude: " << props.magneticFieldAmplitude << std::endl;
			out << "\tCharge density: " << props.rhoInit << std::endl;
			return out;
		}

	};

}
