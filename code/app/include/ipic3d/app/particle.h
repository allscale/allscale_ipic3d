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

		// user-specified default constructor to ensure proper initialization
		Particle() : position(), velocity(), q(), mass(), velocityStar() {};

		void updatePosition(double dt) {
			position += velocity * dt;
		}

		void updateVelocity(const Force& force, double dt) {
			velocity += force * dt / mass;
		}

		void updateVelocityBorisStyle(const Force& E, const Force& B, double dt) {

			auto k = q / mass * 0.5 * dt;

			auto t = k * B;

			auto t_mag2 = allscale::api::user::data::sumOfSquares(t);

			auto s = (2.0 * t) / (1+t_mag2);

			auto v_minus = velocity + k * E;

			auto v_prime = v_minus + crossProduct(v_minus,t);

			auto v_plus = v_minus + crossProduct(v_prime,s);

			velocity = v_plus + k * E;

		}
		
		friend std::ostream& operator<<(std::ostream& out, const Particle& p) {
			out << "Particle: " << std::endl;
			out << "\tPosition: " << p.position << std::endl;
			out << "\tVelocity: " << p.velocity << std::endl;
			out << "\tVelocityStar: " << p.velocityStar << std::endl;
			out << "\tCharge: " << p.q << std::endl;
			out << "\tMass: " << p.mass << std::endl;
			return out;
		}

	};

} // end namespace ipic3d
