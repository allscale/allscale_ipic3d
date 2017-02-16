#pragma once

#include "ipic3d/app/vector.h"

namespace ipic3d {


	using Force = Vector3<double>;

	struct Particle {

		Vector3<double> position;			// position (absolute - TODO: relative to cell center)

		Vector3<double> velocity;			// velocity of this particle

		double q;							// charge of this particle

		double mass;						// the mass of this particle

		Vector3<double> velocityStar;  	// auxiliary parameters for the Boris mover

		void updatePosition(double dt) {
			position += velocity * dt;
		}

		void updateVelocity(const Force& force, double dt) {
			velocity += force * dt / mass;
		}

		void updateVelocityBorisStyle(const Force& E, const Force& B, double dt) {
			Vector3<double> v = velocity;

			auto k = q / mass * 0.5 * dt;

			auto t = k * B;

			auto t_mag2 = sumOfSquares(t);

			auto s = (2.0 * t) / (1+t_mag2);

			auto v_minus = v + k * E;

			auto v_prime = v_minus + crossProduct(v_minus,t);

			auto v_plus = v_minus + crossProduct(v_prime,s);

			v = v_plus + k * E;

			velocity = v;
		}

	};

} // end namespace ipic3d
