#include <iostream>

#include "matrix2D.h"

template<typename T>
void title(const T& x)
{
	std::cout << '\n' << x << "\n--------------------" << std::endl;
}

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

	title("x");
	x.print(2);
	title("|x|");
	cout << y << endl;
	title("x^-1");
	z.print(2);
	title("x * x^-1 = I");
	I.print(2);

	matrix2D<double> b(3,1);
	b[0]={1};
	b[1]={1};
	b[2]={1};

	lu_decomposition<double> LU(x);

	title("LU(x)");
	LU.print(2);
	title("|x|");
	cout << x.det() << endl;
	title("|LU(x)|");
	cout << LU.det() << endl;

	auto solution = LU.solve(b);
	title("b");
	b.print(2);
	title("LU(x) solve for b");
	solution.print(2);

	auto matrixExp = x.expm();
	title("expm(x)");
	matrixExp.print(2);

	solution = x.solveL2(b);
	title("L2 norm solve x for b");
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

	title("A");
	A.print();
	title("b");
	b.print();
	title("L2 norm solve A for b");
	solution = A.solveL2(b);
	solution.print(2);
	cout << endl;

	char _1;
	cin >> _1;

	return 0;
}
