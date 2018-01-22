#include <iostream>

#include "matrix2D.h"

int main(int argc, char** argv)
{
	using namespace linear_algebra;
	using std::cin;
	using std::cout;
	using std::endl;

	matrix2D<double> x(3, 3);

	x[0] = { 3, 0, 2 };
	x[1] = { 2, 0,-2 };
	x[2] = { 0, 1, 1 };

	auto y = x.det();

	auto z = x.inverse();

	auto I = x * z;

	x.print(2);
	cout << endl;
	cout << y << endl;
	cout << endl;
	z.print(2);
	cout << endl;
	I.print(2);
	cout << endl;

	matrix2D<double> b(3,1);
	b[0]={1};
	b[1]={1};
	b[2]={1};
	lu_decomposition<double> LU(x);
	LU.print(2);
	cout << endl;

	cout << x.det() << endl;
	cout << LU.det() << endl;

	matrix2D<double> solution = LU.solve(b);
	solution.print(2);
	cout << endl;

	matrix2D<double> matrixExp = expm(x);
	matrixExp.print(2);
	cout << endl;

	solution = solveL2(x,b);
	solution.print(2);
	cout << endl;

	matrix2D<double> A(4,3);
	A[0] = { 3, 0, 2 };
	A[1] = { 2, 0,-2 };
	A[2] = { 0, 1, 1 };
	A[3] = { 4, 1, 1 };
	b = matrix2D<double>(4,1);
	b[0]={1};
	b[1]={1};
	b[2]={1};
	b[3]={1};

	solution = solveL2(A,b);
	solution.print(2);
	cout << endl;



	char _1;
	cin >> _1;

	return 0;
}
