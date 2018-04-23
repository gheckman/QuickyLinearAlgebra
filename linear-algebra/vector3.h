/**
@file vector3.h
@brief 3D vector.
*/

#pragma once

/**
3D vector.
X is right. Y is forward. Z is up.
@tparam T The type of data each element of the vector holds. Recommended: float or double.
*/
namespace linear_algebra
{
	/**
	3D vector.
	X is right. Y is forward. Z is up.
	@tparam T The type of data each element of the vector holds. Recommended: float or double.
	*/
	template<typename T>
	class vector3
	{
		constexpr static int DIM = 3;

	public:

		/**
		Default constructor.
		Represents a vector3 with components {0, 0, 0}.
		*/
		vector3() : vector3({ 0, 0, 0 }) {}

		/**
		Parameterized constructor.
		Represents a vector3 with components {x, y, z}.
		*/
		vector3(T x, T y, T z) : vector3({ x, y, z }) {}

		/**
		Parameterized constructor.
		@param [in] list Initializer list.
		Represents a vector3 with components {x, y, z} taken from the first 3 elements of the list.
		*/
		vector3(std::initializer_list<T> list) : vector3(list.begin(), list.end()) {}

		/**
		Parameterized constructor.
		Represents a vector3 with components {x, y, z} taken from the first 3 elements of the iterator range.
		@param [in] first Iterator to the first element.
		@param [in] last  Iterator to 1 past the last element.
		@tparam Iter Forward iterator.
		*/
		template <typename Iter>
		vector3(Iter first, Iter last)
		{
			auto dist = std::min(DIM, std::distance(first, last));
			std::copy(first, std::advance(first, dist), _xyz.begin());
		}

		/**
		Adds 2 vector3s.
		@tparam U Some vector3 element type.
		@param [in] rhs vector3 on the right hand side.
		@return The elementwise sum of the 2 vector3s.
		*/
		template<typename U>
		auto operator += (const vector3<U>& rhs)
		{
			*this.x() += rhs.x();
			*this.y() += rhs.y();
			*this.z() += rhs.z();
			return *this;
		}

		/**
		Adds 2 vector3s.
		@tparam U Some vector3 element type.
		@param [in] rhs vector3 on the right hand side.
		@return The elementwise sum of the 2 vector3s.
		*/
		template<typename U>
		auto operator + (const vector3<U>& rhs) const
		{
			auto temp(*this);
			return temp += rhs;
		}

		/**
		Subtracts 2 vector3s.
		@tparam U Some vector3 element type.
		@param [in] rhs vector3 on the right hand side.
		@return The elementwise difference of the 2 vector3s.
		*/
		template<typename U>
		auto operator -= (const vector3<U>& rhs)
		{
			*this.x() -= rhs.x();
			*this.y() -= rhs.y();
			*this.z() -= rhs.z();
			return *this;
		}

		/**
		Subtracts 2 vector3s.
		@tparam U Some vector3 element type.
		@param [in] rhs vector3 on the right hand side.
		@return The elementwise difference of the 2 vector3s.
		*/
		template<typename U>
		auto operator - (const vector3<U>& rhs) const
		{
			auto temp(*this);
			return temp -= rhs;
		}

		/**
		Multiplies a vector3 by a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return The product of the vector3 with the scalar.
		*/
		template<typename U>
		auto& operator *= (const U& rhs)
		{
			_xyz[0] *= rhs;
			_xyz[1] *= rhs;
			_xyz[2] *= rhs;
			return *this;
		}

		/**
		Multiplies a vector3 by a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return The product of the vector3 with the scalar.
		*/
		template<typename U>
		auto operator * (const U& rhs) const
		{
			auto temp(*this);
			return temp *= rhs;
		}

		/**
		Divides a vector3 by a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return The product of the vector3 with the inverse of the scalar.
		*/
		template<typename U>
		auto& operator /= (const U& rhs)
		{
			return *this * (1 / rhs);
		}

		/**
		Divides a vector3 by a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return The product of the vector3 with the inverse of the scalar.
		*/
		template<typename U>
		auto operator / (const U& rhs) const
		{
			auto temp(*this);
			return temp /= rhs;
		}

		/**
		Negation operator.
		@return Elementwise negation of the vector3.
		*/
		auto operator -() const
		{
			return vector3(-x(), -y(), -z());
		}

		/**
		Dot product of the 2 vector3s.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return The product of the vector3 with the inverse of the scalar.
		*/
		template<typename U>
		auto dot(const vector3<U>& rhs)
		{
			return x() * rhs.x() + y() * rhs.y() + z() * rhs.z();
		}

		/**
		Cross product of the 2 vector3s.
		@tparam U Some vector3 element type.
		@param [in] rhs vector3 on the right hand side.
		@return The product of the vector3 with the inverse of the scalar.
		*/
		template<typename U>
		auto cross(const vector3<U>& rhs)
		{
			using ResultT = decltype(std::declval<T>() * std::declval<U>());
			return vector3<ResultT>
				(
					y() * rhs.z() - z() * rhs.y(),
					z() * rhs.x() - x() * rhs.z(),
					x() * rhs.y() - y() * rhs.x()
					);
		}

		/**
		Gets the vector3's L2 norm.
		This is equal to the magnitude squared.
		Prefered over normalized and norm when possible.
		@return The vector3's L2 norm.
		*/
		auto l2norm() const
		{
			return x() * x() + y() * y() + z() * z();
		}

		/**
		Gets the vector3's norm.
		Also known as magnitude.
		@return The quaternion's magnitude.
		*/
		auto norm() const
		{
			return ::sqrt(l2norm());
		}

		/**
		Normalizes the vector3.
		Divides the vector3 by its norm.
		@return The normalized vector3.
		*/
		auto normalized() const
		{
			return *this / magnitude();
		}

		/**
		Outputs a formatted text representation of the vector3.
		@param [in, out] os Output stream to write the formatted text to.
		@param [in]      v  vector3 to be written to the output stream.
		@return The output stream.
		*/
		friend auto& operator << (std::ostream& os, const vector3<T>& v)
		{
			return os << "[" << v.x() << ", " << v.y() << ", " << v.z() << "]";
		}

		/**
		Gets the vector associated with the direction 'right'.
		@return The vector associated with the direction 'right'.
		*/
		auto inline static right() { return vector3(T{ 0 }, T{ 0 }, T{ 1 }); }

		/**
		Gets the vector associated with the direction 'forward'.
		@return The vector associated with the direction 'forward'.
		*/
		auto inline static forward() { return vector3(T{ 0 }, T{ 1 }, T{ 0 }); }

		/**
		Gets the vector associated with the direction 'up'.
		@return The vector associated with the direction 'up'.
		*/
		auto inline static up() { return vector3(T{ 0 }, T{ 0 }, T{ 1 }); }

		/**
		Gets the vector associated with the direction 'left'.
		@return The vector associated with the direction 'left'.
		*/
		auto inline static left() { return -right(); }

		/**
		Gets the vector associated with the direction 'back'.
		@return The vector associated with the direction 'back'.
		*/
		auto inline static back() { return -forward(); }

		/**
		Gets the vector associated with the direction 'down'.
		@return The vector associated with the direction 'down'.
		*/
		auto inline static down() { return -up(); }

		/**
		Sets the X component of the vector3.
		*/
		inline void x(T x) { _xyz[0] = x; }

		/**
		Sets the Y component of the vector3.
		*/
		inline void y(T y) { _xyz[1] = y; }

		/**
		Sets the Z component of the vector3.
		*/
		inline void z(T z) { _xyz[2] = z; }

		/**
		Gets the X component of the vector3.
		@return X component of the vector3.
		*/
		inline T x() const { return _xyz[0]; }

		/**
		Gets the Y component of the vector3.
		@return Y component of the vector3.
		*/
		inline T y() const { return _xyz[1]; }

		/**
		Gets the Z component of the vector3.
		@return Z component of the vector3.
		*/
		inline T z() const { return _xyz[2]; }

	private:
		std::array<T, DIM> _xyz;
	};
}