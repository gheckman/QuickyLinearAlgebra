#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <vector>

namespace linear_algebra
{
	template <typename T>
	class matrix_row
	{
	public:
		matrix_row(std::vector<T>& row) :
			row_(row)
		{}

		T& operator [] (const size_t i)
		{
			return row_[i];
		}

		const T& operator [] (const size_t i) const
		{
			return row_[i];
		}

		template <typename InputIter>
		void assign(InputIter first, InputIter last)
		{
			auto dist = std::min(std::distance(cbegin(row_), cend(row_)), std::distance(first, last));
			auto real_last = first;
			std::advance(real_last, dist);
			std::transform(first, real_last, begin(row_), [](T x) { return x; });
		}

		void assign(const std::initializer_list<T>& list)
		{
			assign(cbegin(list), cend(list));
		}

		matrix_row& operator = (const std::initializer_list<T>& list)
		{
			assign(list);
			return *this;
		}

	private:
		std::vector<T>& row_;
	};

	template <typename T>
	class matrix2D
	{
	public:
		matrix2D(size_t rows, size_t cols) :
			rows_(rows),
			cols_(cols)
		{
			data_.assign(rows, std::vector<T>{});
			std::for_each(begin(data_), end(data_), [=](std::vector<T>& x) { x.assign(cols, 0); });
		}

#pragma region operators

#pragma region access

		matrix_row<T> operator [] (const size_t i)
		{
			return data_.at(i);
		}

		const matrix_row<T> operator [] (const size_t i) const
		{
			// Dangerous to take away const, but I'm adding back in, so...
			auto& non_const_vector = const_cast<std::vector<T>&>(data_.at(i));
			const matrix_row<T> row(non_const_vector);
			return row;
		}

#pragma endregion

#pragma region matrix math

		matrix2D& operator *= (const matrix2D& rhs)
		{
			auto new_matrix(*this);

			*this = new_matrix * rhs;

			return *this;
		}

		matrix2D operator * (const matrix2D& rhs) const
		{
			matrix2D new_matrix(rows_, rhs.cols_);

			assert(cols_ == rhs.rows_);

			for (size_t i = 0; i < rows_; ++i)
				for (size_t j = 0; j < rhs.cols_; ++j)
					for (size_t k = 0; k < cols_; ++k)
						new_matrix.data_[i][j] += data_[i][k] * rhs.data_[k][j];

			return new_matrix;
		}

		matrix2D& operator += (const matrix2D& rhs)
		{
			return add_or_sub(rhs, std::plus<T>());
		}

		matrix2D operator + (const matrix2D& rhs) const
		{
			auto new_matrix(*this);
			new_matrix += rhs;
			return new_matrix;
		}

		matrix2D& operator -= (const matrix2D& rhs)
		{
			return add_or_sub(rhs, std::minus<T>());
		}

		matrix2D operator - (const matrix2D& rhs) const
		{
			auto new_matrix(*this);
			new_matrix -= rhs;
			return new_matrix;
		}

#pragma endregion

#pragma region scalar math

		matrix2D& operator *= (const T& rhs)
		{
			return scalar_op(rhs, std::multiplies<T>());
		}

		matrix2D operator * (const T& rhs) const
		{
			auto new_matrix(*this);
			new_matrix *= rhs;
			return new_matrix;
		}

		matrix2D& operator /= (const T& rhs)
		{
			return scalar_op(rhs, std::divides<T>());
		}

		matrix2D operator / (const T& rhs) const
		{
			auto new_matrix(*this);
			new_matrix /= rhs;
			return new_matrix;
		}

		matrix2D& operator += (const T& rhs)
		{
			return scalar_op(rhs, std::plus<T>());
		}

		matrix2D operator + (const T& rhs) const
		{
			auto new_matrix(*this);
			new_matrix += rhs;
			return new_matrix;
		}

		matrix2D& operator -= (const T& rhs)
		{
			return scalar_op(rhs, std::minus<T>());
		}

		matrix2D operator - (const T& rhs) const
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
}