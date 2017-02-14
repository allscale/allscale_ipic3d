#pragma once

#include "ipic3d/vector.h"

namespace ipic3d {


	using Force = Vector3<double>;

	struct Particle {

		double x,y,z;						// position (absolute - TODO: relative to cell center)

		double dx,dy,dz;					// velocity

		double q;							// charge divided over the mass of spacies

		double mass;						// the mass of this particle

		double vxstar, vystar, vzstar;  	// auxiliary parameters for the Boris mover

		void updatePosition(double dt) {
			x += dx * dt;
			y += dy * dt;
			z += dz * dt;
		}

		void updateVelocity(const Force& force, double dt) {
			dx += force.x * dt / mass;
			dy += force.y * dt / mass;
			dz += force.z * dt / mass;
		}

		void updateVelocityBorisStyle(const Force& E, const Force& B, double dt) {
			Vector3<double> v { dx, dy, dz };

			auto k = q / mass * 0.5 * dt;

			auto t = k * B;

			auto t_mag2 = t.x*t.x + t.y*t.y + t.z*t.z;

			auto s = (2.0 * t) / (1+t_mag2);

			auto v_minus = v + k * E;

			auto v_prime = v_minus + cross(v_minus,t);

			auto v_plus = v_minus + cross(v_prime,s);

			v = v_plus + k * E;

			dx = v.x;
			dy = v.y;
			dz = v.z;
		}

	};


	struct Particle2 {

		Vector3<double> x;				// the position

		Vector3<double> v;				// the velocity

		double q;						// the charge

		double m;						// the mass


		void updatePosition(double dt) {
			x += v * dt;
		}

		void updateVelocity(const Force& force, double dt) {
			v += force * dt / m;
		}

		Force getForce(const Vector3<double>& E, const Vector3<double>& B) const {
			return q * (E + cross(v,B));
		}

	};


} // end namespace ipic3d
