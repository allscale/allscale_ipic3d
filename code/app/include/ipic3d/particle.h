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
