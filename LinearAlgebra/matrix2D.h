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

		size_t rowSize() const {return rows_;}
		size_t columnSize() const {return cols_;}


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

	template <typename T>
	class lu_decomposition
	{
	public:
		lu_decomposition(matrix2D<T> input):
			parity(1),
			index(input.rowSize(),0),
			tiny_value(1e-20),
			A(input)
		{
			assert(A.rowSize() == A.columnSize());
			decompose();
		}

		void decompose()
		{
			//transcribed from c numerical recipes, as always
			//this destroys the matrix, replacing it with the LU decomposition

			std::vector<T> scaling_vector;

			for(int i = 0; i < A.rowSize(); i++)
			{
				T big = 0.0;
				for (int j = 0; j< A.rowSize(); j++)
				{
					if(abs(A[i][j])>big)
						big = abs(A[i][j]);
				}
				assert(big > 0.0); // singular
				scaling_vector.push_back(1.0/big);

			}

			//crout's method
			int imax = 0;
			for(int j=0; j<A.rowSize(); j++)
			{
				for(int i = 0; i<j; i++)
				{
					T sum=A[i][j];	
					for(int k=0; k<i;k++)
						sum-= A[i][k]*A[k][j];
					A[i][j]=sum;
				}
				//search for pivot
				T biggest_value =0;
				for(int i = j; i < A.rowSize(); i++)
				{
					T sum = A[i][j];
					for (int k=0;k<j;k++)
						sum -= A[i][k]*A[k][j];
					A[i][j]=sum;
					if ( scaling_vector[i]*abs(sum) > biggest_value)
					{
						biggest_value = scaling_vector[i]*(abs(sum));
						imax = i;
					}
				}
				//do pivot
				if(j != imax)
				{
					for(int k=0; k<A.rowSize(); k++)
						std::swap(A[imax][k],A[j][k]);
					parity = -parity;
					std::swap(scaling_vector[j],scaling_vector[imax]);
				}
				index[j]=imax; //keep track of which rows were pivoted
				if(A[j][j]==0.0)
					A[j][j]=tiny_value;  //the attempt has failed, singular matrix.  using tiny_value will prevent inf values, which can salvage some info for postmortem
				for(int i=j+1; i<A.rowSize(); i++)
					A[i][j]/=A[j][j];
			}
		}

		matrix2D<T> solve(const matrix2D<T>& b)
		{
			//make sure the dimensions mesh
			assert(A.columnSize() == b.rowSize());

			matrix2D<T> x = b;
			//solve all columns
			for(int column_i = 0; column_i < x.columnSize(); column_i++)
				{
				int first_nonzero_element = -1;
				for(int i=0; i<A.rowSize(); i++) //forward substitute
				{
					T sum = x[index[i]][column_i];
					if(first_nonzero_element>-1)
					{
						for (int j = first_nonzero_element; j<i;j++)
							sum-=A[i][j]*x[j][column_i];
					}
					else
					{
						if(sum>0.0)
							first_nonzero_element=i;
					}
					x[i][column_i]=sum;
	
				}
				for(int i=A.rowSize()-1; i>=0; i--) //back substitute
				{
					T sum=x[i][column_i];
					for(int j=i+1;j<A.rowSize(); j++)
						sum-=A[i][j]*x[j][column_i];
					x[i][column_i] = sum/A[i][i];
				}
			}
			return x;
		}

		T det()
		{
			T d = parity;
			for(int i=0; i < A.rowSize(); i++)
				d*=A[i][i];
			return d;
		}

		void print() const
		{
			A.print();
		}
		void print(std::streamsize precision) const
		{
			A.print(precision);
		}


	private:
		T tiny_value;
		matrix2D<T> A;
		std::vector<int> index;
		T parity;
	};

	matrix2D<double> eye(size_t matrix_size)
	{
		matrix2D<double> A(matrix_size,matrix_size);
		for(int i=0; i< matrix_size; i++)
		{
			A[i][i]=1.0;
		}
		return A;
	}

	template<typename T>
	matrix2D<T> expm(matrix2D<T> A)
	{
		assert(A.rowSize() == A.columnSize());
		T error = 1.0;
		int shots = 1;

		matrix2D<T> Aexpm=eye(A.rowSize());
		matrix2D<T> A_collected=eye(A.rowSize());
		T determinant_old = 1.0;
		T determinant;

		while(error > 1e-8)
		{
			A_collected=A_collected*A/((T)shots);
			Aexpm+=A_collected;	
			determinant_old=determinant;
			determinant = Aexpm.det();
			error = abs(determinant - determinant_old);
			shots++;
		}
		return Aexpm;
	}

	template<typename T>
	matrix2D<T> solveL2(matrix2D<T> A,matrix2D<T> b)
	{
		//normal equations
		lu_decomposition<T> LU(A.transpose()*A);
		return LU.solve(A.transpose()*b);
	}

}
