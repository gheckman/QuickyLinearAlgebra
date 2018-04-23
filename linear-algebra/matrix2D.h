/**
@file matrix2D.h
@brief It's a matrix.
*/

#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <vector>
#include <cmath>

namespace linear_algebra
{
	/**
	A single row from the matrix.
	Limited API into the storage structure of a single matrix row.
	@tparam T The type of data each element of the matrix row holds. Recommended: float or double.
	*/
	template <typename T>
	class matrix_row
	{
	public:

		/**
		Default constructor.
		Represents a single matrix row.
		*/
		matrix_row(std::vector<T>& row) :
			row_(row)
		{}

		/**
		Element access operator.
		@param [in] i Index into the row. 
		@return Element at the given index.
		*/
		auto& operator [] (const size_t i)
		{
			return row_[i];
		}

		/**
		Element access operator.
		@param [in] i Index into the row.
		@return Element at the given index.
		*/
		const T& operator [] (const size_t i) const
		{
			return row_[i];
		}

		/**
		Assigns an iterator range to the row.
		@tparam FwdIt Forward iterator.
		@param [in] first Iterator to the first element.
		@param [in] last  Iterator to 1 past the last element.
		*/
		template <typename FwdIt>
		void assign(FwdIt first, FwdIt last)
		{
			auto dist = std::min(std::distance(cbegin(row_), cend(row_)), std::distance(first, last));
			auto real_last = first;
			std::advance(real_last, dist);
			// todo [9] std::copy here instead?
			std::transform(first, real_last, begin(row_), [](T x) { return x; });
		}

		/**
		Assigns an initializer list to the row.
		@param [in] list Initializer list.
		*/
		void assign(const std::initializer_list<T>& list)
		{
			assign(cbegin(list), cend(list));
		}

		/**
		Assignment operator.
		Assigns an initializer list to the row.
		@param [in] list Initializer list.
		@return this.
		*/
		auto& operator = (const std::initializer_list<T>& list)
		{
			assign(list);
			return *this;
		}

	private:
		std::vector<T>& row_;
	};

	/**
	Default constructor.
	Zero matrix.
	@param [in] rows Number of rows.
	@param [in] cols Number of columns.
	*/
	template <typename T>
	class matrix2D
	{
	public:
		matrix2D(size_t rows, size_t cols) :
			rows_(rows),
			cols_(cols)
		{
			data_.assign(rows, std::vector<T>{});
			for (auto& row : data_)
				row.assign(cols, 0);
		}

#pragma region operators

#pragma region access

		/**
		Index operator.
		@param [in] i Index into the rows.
		@return Row of the matrix at the given index.
		*/
		auto operator [] (const size_t i)
		{
			return matrix_row<T>(data_.at(i));
		}

		/**
		Index operator.
		@param [in] i Index into the rows.
		@return Row of the matrix at the given index.
		*/
		const auto operator [] (const size_t i) const
		{
			// Dangerous to take away const, but I'm adding back in, so...
			auto& non_const_vector = const_cast<std::vector<T>&>(data_.at(i));
			const matrix_row<T> row(non_const_vector);
			return row;
		}

#pragma endregion

#pragma region matrix math

		/**
		Multiplies 2 matrix2Ds.
		@tparam U Some matrix2D element type.
		@return The product of the matrix2Ds.
		*/
		template<typename U>
		auto& operator *= (const matrix2D<U>& rhs)
		{
			*this = *this * rhs;
			return *this;
		}

		/**
		Multiplies 2 matrix2Ds.
		@tparam U Some matrix2D element type.
		@return The product of the matrix2Ds.
		*/
		template<typename U>
		auto operator * (const matrix2D<U>& rhs) const
		{
			using ResultT = decltype(std::declval<T>() * std::declval<U>());
			matrix2D<ResultT> new_matrix(rows_, rhs.cols_);

			assert(cols_ == rhs.rows_);

			for (size_t i = 0; i < rows_; ++i)
				for (size_t j = 0; j < rhs.cols_; ++j)
					for (size_t k = 0; k < cols_; ++k)
						new_matrix.data_[i][j] += data_[i][k] * rhs.data_[k][j];

			return new_matrix;
		}

		/**
		Adds 2 matrix2Ds.
		@tparam U Some matrix2D element type.
		@return Elementwise addition of the matrix2Ds.
		*/
		template<typename U>
		auto& operator += (const matrix2D<U>& rhs)
		{
			return add_or_sub(rhs, std::plus<T>());
		}

		/**
		Adds 2 matrix2Ds.
		@tparam U Some matrix2D element type.
		@return Elementwise addition of the matrix2Ds.
		*/
		template<typename U>
		auto operator + (const matrix2D<U>& rhs) const
		{
			auto new_matrix(*this);
			new_matrix += rhs;
			return new_matrix;
		}

		/**
		Subtracts 2 matrix2Ds.
		@tparam U Some matrix2D element type.
		@return Elementwise subtraction of the matrix2Ds.
		*/
		template<typename U>
		auto& operator -= (const matrix2D<U>& rhs)
		{
			return add_or_sub(rhs, std::minus<T>());
		}

		/**
		Subtracts 2 matrix2Ds.
		@tparam U Some matrix2D element type.
		@return Elementwise subtraction of the matrix2Ds.
		*/
		template<typename U>
		auto operator - (const matrix2D<U>& rhs) const
		{
			auto new_matrix(*this);
			new_matrix -= rhs;
			return new_matrix;
		}

#pragma endregion

#pragma region scalar math

		/**
		Multiplies a matrix2D by a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return The product of the matrix2D with the scalar.
		*/
		template<typename U>
		auto& operator *= (const U& rhs)
		{
			return scalar_op(rhs, std::multiplies<U>());
		}

		/**
		Multiplies a matrix2D by a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return The product of the matrix2D with the scalar.
		*/
		auto operator * (const T& rhs) const
		{
			auto new_matrix(*this);
			new_matrix *= rhs;
			return new_matrix;
		}

		/**
		Divides a matrix2D by a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return The product of the matrix2D with the inverse of the scalar.
		*/
		auto& operator /= (const T& rhs)
		{
			return scalar_op(rhs, std::divides<T>());
		}

		/**
		Divides a matrix2D by a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return The product of the matrix2D with the inverse of the scalar.
		*/
		auto operator / (const T& rhs) const
		{
			auto new_matrix(*this);
			new_matrix /= rhs;
			return new_matrix;
		}

		/**
		Adds a matrix2D and a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return Addition of each matrix element with the scalar.
		*/
		auto& operator += (const T& rhs)
		{
			return scalar_op(rhs, std::plus<T>());
		}

		/**
		Adds a matrix2D and a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return Addition of each matrix element with the scalar.
		*/
		auto operator + (const T& rhs) const
		{
			auto new_matrix(*this);
			new_matrix += rhs;
			return new_matrix;
		}

		/**
		Subtracts a matrix2D and a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return Addition of each matrix element with the negation of the scalar.
		*/
		auto& operator -= (const T& rhs)
		{
			return scalar_op(rhs, std::minus<T>());
		}

		/**
		Subtracts a matrix2D and a scalar.
		@tparam U Some scalar type.
		@param [in] rhs Scalar on the right hand side.
		@return Addition of each matrix element with the negation of the scalar.
		*/
		auto operator - (const T& rhs) const
		{
			auto new_matrix(*this);
			new_matrix -= rhs;
			return new_matrix;
		}

#pragma endregion

#pragma endregion

#pragma region math non-operator

		matrix2D transpose() const
		{
			matrix2D new_matrix(cols_, rows_);

			for (size_t i = 0; i < rows_; ++i)
				for (size_t j = 0; j < cols_; ++j)
					new_matrix.data_[j][i] = data_[i][j];

			return new_matrix;
		}

		T det()
		{
			T d{};

			assert(rows_ == cols_);

			auto n = rows_;

			if (n == 2)
				return data_[0][0] * data_[1][1] - data_[1][0] * data_[0][1];
			if (n == 1)
				return data_[0][0];

			T sign{ 1 };

			for (size_t i = 0; i < n; i++)
			{
				auto temp = cofactor(0, i);
				d += sign * data_[0][i] * temp.det();
				sign = -sign;
			}

			return d;
		}

		matrix2D cofactor(int p, int q)
		{
			matrix2D new_matrix(rows_ - 1, cols_ - 1);

			assert(rows_ == cols_);

			auto n = rows_;

			int i = 0, j = 0;

			for (size_t row = 0; row < n; ++row)
			{
				for (size_t col = 0; col < n; ++col)
				{
					if (row != p && col != q)
					{
						new_matrix.data_[i][j] = data_[row][col];
						j++;

						if (j == n - 1)
						{
							j = 0;
							i++;
						}
					}
				}
			}

			return new_matrix;
		}

		matrix2D minor()
		{
			matrix2D new_matrix(rows_, cols_);

			assert(rows_ == cols_);

			for (size_t i = 0; i < rows_; ++i)
				for (size_t j = 0; j < rows_; ++j)
					new_matrix[i][j] = cofactor(i, j).det();

			return new_matrix;
		}

		matrix2D cofactor_matrix()
		{
			return minor().checkerboard_negate();
		}

		matrix2D checkerboard_negate()
		{
			matrix2D new_matrix(*this);

			assert(rows_ == cols_);

			T sign{ 1 };

			for (size_t i = 0; i < rows_; ++i)
			{
				for (size_t j = 0; j < rows_; ++j)
				{
					new_matrix[i][j] = sign * new_matrix[i][j];
					sign = -sign;
				}
				if (rows_ % 2 == 0)
					sign = -sign;
			}

			return new_matrix;
		}

		matrix2D inverse()
		{
			assert(rows_ == cols_);

			matrix2D cof(cofactor_matrix());
			T det{};

			// fast determinate calculation since I needed the cofactor matrix anyway
			for (size_t i = 0; i < rows_; ++i)
				det += data_[0][i] * cof.data_[0][i];

			return cof.transpose() / det;
		}

		matrix2D<T> expm()
		{
			// NRVO delaration
			auto expm_matrix = identity(row_size());

			constexpr double error_expm = 1e-8;

			assert(row_size() == column_size());

			T error{ 1.0 };
			int shots{ 1 };

			auto collected = identity(row_size());
			T determinant_old{ 1.0 };
			T determinant{};

			while (error > error_expm)
			{
				collected *= *this / (T)shots;
				expm_matrix += collected;
				determinant_old = determinant;
				determinant = expm_matrix.det();
				error = std::abs(determinant - determinant_old);
				++shots;
			}

			return expm_matrix;
		}

		matrix2D<T> solveL2(const matrix2D<T>& b)
		{
			//normal equations
			lu_decomposition<T> LU(transpose() * *this);
			return LU.solve(transpose() * b);
		}

#pragma endregion

		void print(std::streamsize precision) const
		{
			std::cout << std::fixed << std::setprecision(precision);

			for (const auto& row : data_)
			{
				for (const auto& t : row)
				{
					std::cout << t << ' ';
				}
				std::cout << std::endl;
			}
		}

		void print() const
		{
			print(6);
		}

		size_t row_size() const { return rows_; }
		size_t column_size() const { return cols_; }


	private:
		size_t rows_;
		size_t cols_;
		std::vector<std::vector<T>> data_;

		template <typename BinaryOperation>
		matrix2D& add_or_sub(const matrix2D& rhs, const BinaryOperation& op)
		{
			assert(cols_ == rhs.cols_ && rows_ == rhs.rows_);

			for (size_t i = 0; i < rows_; ++i)
				for (size_t j = 0; j < rhs.cols_; ++j)
					data_[i][j] = op(data_[i][j], rhs.data_[i][j]);

			return *this;
		}

		template <typename BinaryOperation>
		matrix2D& scalar_op(const T& rhs, const BinaryOperation& op)
		{
			for (size_t i = 0; i < rows_; ++i)
				for (size_t j = 0; j < cols_; ++j)
					data_[i][j] = op(data_[i][j], rhs);

			return *this;
		}
	};

	matrix2D<double> identity(size_t matrix_size)
	{
		matrix2D<double> I(matrix_size, matrix_size);
		for (size_t i = 0; i < matrix_size; ++i)
			I[i][i] = 1.0;
		return I;
	}

	template <typename T>
	class lu_decomposition
	{
	public:
		lu_decomposition(matrix2D<T> input) :
			parity_(1),
			index_(input.row_size(), 0),
			tiny_value_(1e-20),
			A_(input)
		{
			assert(A_.row_size() == A_.column_size());
			decompose();
		}

		matrix2D<T> solve(const matrix2D<T>& b)
		{
			//make sure the dimensions mesh
			assert(A_.column_size() == b.row_size());

			auto x = b;
			//solve all columns
			for (size_t column_i = 0; column_i < x.column_size(); ++column_i)
			{
				int first_nonzero_element = -1;
				for (size_t i = 0; i < A_.row_size(); ++i) //forward substitute
				{
					T sum = x[index_[i]][column_i];
					if (first_nonzero_element > -1)
					{
						for (size_t j = first_nonzero_element; j < i; ++j)
							sum -= A_[i][j] * x[j][column_i];
					}
					else
					{
						if (sum > 0.0)
							first_nonzero_element = i;
					}
					x[i][column_i] = sum;

				}
				for (int i = A_.row_size() - 1; i >= 0; --i) //back substitute
				{
					T sum = x[i][column_i];
					for (size_t j = i + 1; j < A_.row_size(); ++j)
						sum -= A_[i][j] * x[j][column_i];
					x[i][column_i] = sum / A_[i][i];
				}
			}
			return x;
		}

		T det()
		{
			T d = parity_;
			for (size_t i = 0; i < A_.row_size(); i++)
				d *= A_[i][i];
			return d;
		}

		void print() const
		{
			A_.print();
		}

		void print(std::streamsize precision) const
		{
			A_.print(precision);
		}

	private:
		T tiny_value_;
		matrix2D<T> A_;
		std::vector<int> index_;
		T parity_;

		void decompose()
		{
			//transcribed from c numerical recipes, as always
			//this destroys the matrix, replacing it with the LU decomposition

			std::vector<T> scaling_vector;

			for (size_t i = 0; i < A_.row_size(); i++)
			{
				T big = 0.0;
				for (size_t j = 0; j < A_.row_size(); j++)
				{
					if (abs(A_[i][j]) > big)
						big = abs(A_[i][j]);
				}
				assert(big > 0.0); // singular
				scaling_vector.push_back(1.0 / big);

			}

			//crout's method
			int imax = 0;
			for (size_t j = 0; j < A_.row_size(); j++)
			{
				for (size_t i = 0; i < j; i++)
				{
					T sum = A_[i][j];
					for (size_t k = 0; k < i; k++)
						sum -= A_[i][k] * A_[k][j];
					A_[i][j] = sum;
				}
				//search for pivot
				T biggest_value = 0;
				for (size_t i = j; i < A_.row_size(); i++)
				{
					T sum = A_[i][j];
					for (size_t k = 0; k < j; k++)
						sum -= A_[i][k] * A_[k][j];
					A_[i][j] = sum;
					if (scaling_vector[i] * abs(sum) > biggest_value)
					{
						biggest_value = scaling_vector[i] * abs(sum);
						imax = i;
					}
				}
				//do pivot
				if (j != imax)
				{
					for (size_t k = 0; k < A_.row_size(); k++)
						std::swap(A_[imax][k], A_[j][k]);
					parity_ = -parity_;
					std::swap(scaling_vector[j], scaling_vector[imax]);
				}
				index_[j] = imax; //keep track of which rows were pivoted
				if (A_[j][j] == 0.0)
					A_[j][j] = tiny_value_;  //the attempt has failed, singular matrix.  using tiny_value will prevent inf values, which can salvage some info for postmortem
				for (size_t i = j + 1; i < A_.row_size(); i++)
					A_[i][j] /= A_[j][j];
			}
		}
	};

}
