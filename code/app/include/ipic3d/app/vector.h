#pragma once

#include "allscale/api/user/data/vector.h"

namespace ipic3d {

	//template<typename T, int dims>
	//using Vector = allscale::api::user::data::Vector<T, dims>;

	//template<typename T>
	//class Vector3 : public Vector<T,3> {
	//public:
	//	T& x = this->at(0);
	//	T& y = this->at(1);
	//	T& z = this->at(2);

	//	Vector3(T x, T y, T z) : x(x), y(y), z(z) {}

	//	Vector3() : x(T()), y(T()), z(T()) {}

	//	Vector3& operator=(const Vector3<T>& other) {
	//		x = other.x;
	//		y = other.y;
	//		z = other.z;
	//		return *this;
	//	}

	//	Vector3& operator+=(const Vector3<T>& other) {
	//		x += other.x;
	//		y += other.y;
	//		z += other.z;
	//		return *this;
	//	}

	//	Vector3& operator-=(const Vector3<T>& other) {
	//		x -= other.x;
	//		y -= other.y;
	//		z -= other.z;
	//		return *this;
	//	}
	//};

	template<typename T>
	struct Vector3 {
		T x;
		T y;
		T z;

		Vector3& operator+=(const Vector3<T>& other) {
			x += other.x;
			y += other.y;
			z += other.z;
			return *this;
		}

		Vector3& operator-=(const Vector3<T>& other) {
			x -= other.x;
			y -= other.y;
			z -= other.z;
			return *this;
		}

	};

	template<typename T>
	Vector3<T> operator+(const Vector3<T>& a, const Vector3<T>& b) {
		return {
			a.x + b.x,
			a.y + b.y,
			a.z + b.z
		};
	}

	template<typename T>
	Vector3<T> operator-(const Vector3<T>& a, const Vector3<T>& b) {
		return{
			a.x - b.x,
			a.y - b.y,
			a.z - b.z
		};
	}


	template<typename T>
	Vector3<T> operator*(Vector3<T> vec, const T& factor) {
		vec.x = vec.x * factor;
		vec.y = vec.y * factor;
		vec.z = vec.z * factor;
		return vec;
	}

	template<typename T>
	Vector3<T> operator*(const T& factor, Vector3<T> vec) {
		vec.x = factor * vec.x;
		vec.y = factor * vec.y;
		vec.z = factor * vec.z;
		return vec;
	}

	template<typename T>
	Vector3<T> operator/(Vector3<T> vec, const T& divisor) {
		vec.x = vec.x / divisor;
		vec.y = vec.y / divisor;
		vec.z = vec.z / divisor;
		return vec;
	}

	template<typename T>
	Vector3<T> crossProduct(const Vector3<T>& a, const Vector3<T>& b) {
		return Vector3<T>{
			a.y * b.z - a.z * b.y,
			a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x
		};
	}

    template <typename T>
    Vector3<T> entrywiseProduct(const Vector3<T>& a, const Vector3<T>& b) {
	    return Vector3<T>{a.x * b.x, a.y * b.y, a.z * b.z};
    }

    template <typename T>
    T sumOfSquares(const Vector3<T>& a) {
	    return a.x * a.x + a.y * a.y + a.z * a.z;
    }


} // end namespace ipic3d
