#pragma once

#include "ipic3d/app/vector.h"

#include <allscale/utils/serializer.h>

namespace ipic3d {

	using Force = Vector3<double>;

	struct Particle : public allscale::utils::trivially_serializable {

		Vector3<double> position;			// position (absolute - TODO: relative to cell center)

		Vector3<double> velocity;			// velocity of this particle

		double q;							// charge of this particle
		double qom;							// charge over mass for this particle

		// user-specified default constructor to ensure proper initialization
		Particle() : position(), velocity(), q(), qom() {};

		/**
 		 * Update position
 		 */
		void updatePosition(double dt) {
        	position += velocity * dt;
        }

		/**
 		 * Update velocity using the Boris method
 		 */ 	
		void updateVelocity(const Force& E, const Force& B, double dt) {

			auto k = qom * 0.5 * dt;

			auto t = k * B;

			auto t_mag2 = allscale::utils::sumOfSquares(t);

			auto s = (2.0 * t) / (1.0 + t_mag2);

			auto v_minus = velocity + k * E;

			auto v_prime = v_minus + crossProduct(v_minus,t);

			auto v_plus = v_minus + crossProduct(v_prime,s);

			velocity = v_plus + k * E;

		}

		friend std::ostream& operator<<(std::ostream& out, const Particle& p) {
			out << "Particle: " << std::endl;
			out << "\tPosition: " << p.position << std::endl;
			out << "\tVelocity: " << p.velocity << std::endl;
			out << "\tCharge: " << p.q << std::endl;
			out << "\tCharge over mass: " << p.qom << std::endl;
			return out;
		}

	};

} // end namespace ipic3d

