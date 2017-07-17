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

		// magnetic field amplitude
		Vector3<double> magneticFieldAmplitude;

		InitProperties(const unsigned numSteps = 1, const std::vector<Vector3<unsigned>>& particlesPerCell = {},
			const std::vector<Vector3<double>>& driftVelocity = {}, const Vector3<double>& magneticFieldAmplitude = { 0,0,0 })
			: numSteps(numSteps), particlesPerCell(particlesPerCell), driftVelocity(driftVelocity), magneticFieldAmplitude(magneticFieldAmplitude) {}

		InitProperties(const Parameters& params) {

			numSteps = params.ncycles;

			for(int i = 0; i < params.ns; i++) {
				driftVelocity.push_back({ params.u0.at(i), params.v0.at(i), params.w0.at(i) });
			}

			for(int i = 0; i < (params.ns); i++) {
				particlesPerCell.push_back({ params.npcelx.at(i), params.npcely.at(i), params.npcelz.at(i) });
			}

			magneticFieldAmplitude = { params.B0.x, params.B0.y, params.B0.z };
		}

	    friend std::ostream& operator<<(std::ostream& out, const InitProperties& props) {
			out << "InitProperties:" << std::endl;
			out << "\tNumber of time steps: " << props.numSteps << std::endl;
			out << "\tNumber of particles per cell: " << props.particlesPerCell << std::endl;
			out << "\tDrift velocity: " << props.driftVelocity << std::endl;
			out << "\tMagnetic field amplitude: " << props.magneticFieldAmplitude << std::endl;
			return out;
		}

	};

}
