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

	char _1;
	cin >> _1;

	return 0;
}