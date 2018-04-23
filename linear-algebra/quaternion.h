/**
	@file quaternion.h
	@brief It's a quaternion.
*/

#pragma once

#include <array>

#include "vector3.h"

/**
It's a quaternion.
@tparam T The type of data each element of the quaternion holds. Recommended: float or double.
*/
namespace linear_algebra
{

	template<typename T>
	class quaternion
	{
		constexpr static int DIM = 4;

	public:

		/**
		Default constructor.
		Represents a normalized quaternion with no rotation.
		*/
		quaternion() : quaternion({ 1, 0, 0, 0 }) {}

		/**
		Parameterized constructor.
		Represents a quaternion with components {w, x, y, z}.
		The quaternion is normalized.
		*/
		quaternion(T w, T x, T y, T z) : quaternion({ w, x, y, z }) {}

		/**
		Parameterized constructor.
		Represents a quaternion with components {w, x, y, z} taken from the first 4 elements of the list.
		The quaternion is normalized.
		*/
		quaternion(std::initializer_list<T> list) : quaternion(list.begin(), list.end()) {}

		/**
		Parameterized constructor.
		Represents a quaternion with components {w, x, y, z} taken from the first 4 elements of the iterator range.
		The quaternion is normalized.
		@tparam Iter Forward iterator.
		*/
		template <typename Iter>
		quaternion(Iter begin, Iter end)
		{
			auto dist = std::min(DIM, std::distance(begin, end));
			std::copy(begin, begin + dist, _wxyz.begin());
			*this = normalized();
		}

		/**
		Multiplies 2 quaternions.
		@tparam U Some quaternion element type.
		@return The product of the 2 quaternions.
		*/
		template<typename U>
		auto operator *= (const quaternion<U>& rhs)
		{
			*this = *this * rhs;
			return *this;
		}

		/**
		Multiplies 2 quaternions.
		@tparam U Some quaternion element type.
		@param [in] rhs quaternion on the right hand side.
		@return The product of the 2 quaternions.
		*/
		template<typename U>
		auto operator * (const quaternion<U>& rhs) const
		{
			using ResultT = decltype(std::declval<T>() * std::declval<U>());
			return quaternion<ResultT>
				(
					w() * rhs.w() - x() * rhs.x() - y() * rhs.y() - z() * rhs.z(),
					w() * rhs.x() + x() * rhs.w() + y() * rhs.z() - z() * rhs.y(),
					w() * rhs.y() - x() * rhs.z() + y() * rhs.w() + z() * rhs.x(),
					w() * rhs.z() + x() * rhs.y() - y() * rhs.x() + z() * rhs.w()
					);
		}

		/**
		Divides 2 quaternions.
		@tparam U Some quaternion element type.
		@return The product of the first quaternion with the conjugate of the second.
		*/
		template<typename U>
		auto operator /= (const quaternion<U>& rhs)
		{
			return *this *= rhs.conjugate() / rhs.l2norm();
		}

		/**
		Divides 2 quaternions.
		@tparam U Some quaternion element type.
		@return The product of the first quaternion with the conjugate of the second.
		*/
		template<typename U>
		auto operator / (const quaternion<U>& rhs) const
		{
			auto temp(*this);
			return temp /= rhs;
		}

		/**
		Multiplies a quaternion by a vector.
		This applies the quaternion's rotation to the vector.
		@tparam U Some vector element type.
		@return The rotated vector.
		*/
		template<typename U>
		auto operator * (const vector3<U>& rhs) const
		{
			using ResultT = typename decltype(T{} *U{});
			quaternion<ResultT> q(0, rhs.x(), rhs.y(), rhs.z());
			q = *this * q * conjugate();
			q *= rhs.norm();
			return vector3<ResultT>(q.x(), q.y(), q.z());
		}

		/**
		Gets the quaternion's conjugate.
		@return The quaternion's conjugate.
		*/
		auto conjugate() const
		{
			return quaternion(w(), -x(), -y(), -z());
		}

		/**
		Gets the quaternion's L2 norm.
		This is equal to the magnitude squared.
		Prefered over normalized and norm when possible.
		@return The quaternion's L2 norm.
		*/
		auto l2norm() const
		{
			return w() * w() + x() * x() + y() * y() + z() * z();
		}

		/**
		Gets the quaternion's norm.
		Also known as magnitude.
		@return The quaternion's magnitude.
		*/
		auto norm() const
		{
			return ::sqrt(l2norm());
		}

		/**
		Normalizes the quaternion.
		Divides the quaternion by its norm.
		@return The normalized quaternion.
		*/
		auto normalized() const
		{
			return *this / norm();
		}

		/**
		Gets a quaternion from the Euler definiton of 3d rotation.
		@param [in] rot Rotation in radians.
		@param [in] v   Axis of rotation.
		@return A quaternion representing the Euler rotation.
		*/
		auto static from_euler(T rot, const vector3<T>& v)
		{
			return from_euler(rot, v.x(), v.y(), v.z());
		}

		/**
		Gets a quaternion from the Euler definiton of 3d rotation.
		@param [in] rot Rotation in radians.
		@param [in] x   Axis of rotation x component.
		@param [in] y   Axis of rotation y component.
		@param [in] z   Axis of rotation z component.
		@return A quaternion representing the Euler rotation.
		*/
		auto static from_euler(T rot, T x, T y, T z)
		{
			auto coef = ::sin(rot / T{ 2 }) / ::sqrt(x * x + y * y + z * z);
			return quaternion(::cos(rot / T{ 2 }), x * coef, y * coef, z * coef);
		}

		/**
		Gets the Euler representation of the quaternion.
		@return The Euler representation of the quaternion.
		*/
		auto to_euler() const
		{
			std::vector<T> rv;
			if (w() == T{ 1 })
				rv = std::vector<T>{ T{ 1 }, T{}, T{}, T{} };
			else
			{
				auto angle = ::atan2(::sqrt(x() * x() + y() * y() + z() * z()), w());
				auto sr = ::sin(angle);
				rv = std::vector<T>{ ::acos(w()) * T { 2 }, x() / sr, y() / sr, z() / sr };
			}
			return rv;
		}

		/**
		Outputs a formatted text representation of the quaternion.
		@param [in, out] os Output stream to write the formatted text to.
		@param [in]      q  quaternion to be written to the output stream.
		@return The output stream.
		*/
		friend auto& operator << (std::ostream& os, const quaternion<T>& q)
		{
			return os << "[" << q.w() << ", " << q.x() << ", " << q.y() << ", " << q.z() << "]";
		}

		/**
		W component of quaternion.
		@return W component of quaternion.
		*/
		inline T w() const { return _wxyz[0]; }

		/**
		X component of quaternion.
		@return X component of quaternion.
		*/
		inline T x() const { return _wxyz[1]; }

		/**
		Y component of quaternion.
		@return Y component of quaternion.
		*/
		inline T y() const { return _wxyz[2]; }

		/**
		Z component of quaternion.
		@return Z component of quaternion.
		*/
		inline T z() const { return _wxyz[3]; }

	private:
		std::array<T, DIM> _wxyz;

		/**
		Multiplies a quaternion by a scalar.
		@tparam U Some scalar type.
		@return The product of the quaternion with the scalar.
		*/
		template<typename U>
		auto& operator *= (const U& rhs)
		{
			_wxyz[0] *= rhs;
			_wxyz[1] *= rhs;
			_wxyz[2] *= rhs;
			_wxyz[3] *= rhs;
			return *this;
		}

		/**
		Multiplies a quaternion by a scalar.
		@tparam U Some scalar type.
		@return The product of the quaternion with the scalar.
		*/
		template<typename U>
		auto operator * (const U& rhs) const
		{
			auto temp(*this);
			return temp *= rhs;
		}

		/**
		Divides a quaternion by a scalar.
		@tparam U Some scalar type.
		@return The product of the quaternion with the inverse of the scalar.
		*/
		template<typename U>
		auto& operator /= (const U& rhs)
		{
			return *this *= (1 / rhs);
		}

		/**
		Divides a quaternion by a scalar.
		@tparam U Some scalar type.
		@return The product of the quaternion with the inverse of the scalar.
		*/
		template<typename U>
		auto operator / (const U& rhs) const
		{
			auto temp(*this);
			return temp /= rhs;
		}
	};
}