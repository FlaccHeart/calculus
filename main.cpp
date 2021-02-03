#include <iostream>
#include "calculus.h"

using namespace std;

int main()
{
	//////////////////////////////
	/// MATRICES
	//////////////////////////////

	// matrices constructor examples
	mat A(3, 3), B(3, 2), C(3, 2), D(3, 2), E(3, 3), F(3, 3), G(3, 3), H(3, 3);

	// random matrix
	A(0, 0) = 1.63; A(0, 1) = 5.73; A(0, 2) = 2.26; 
	A(1, 0) = 2.61; A(1, 1) = 3.86; A(1, 2) = 1.57;
	A(2, 0) = 4.19; A(2, 1) = 6.38; A(2, 2) = 7.28;

	// singular matrix
	B(0, 0) = 12.3; B(0, 1) = 24.3;
	B(1, 0) = 32.1; B(1, 1) = 32.1;
	B(2, 0) = 42.6; B(2, 1) = 15.7;

	// ranks and determinants
	cout << "Ranks and determinants:\n" << endl;
	cout << "rk(A) = " << rk(A) << endl;
	cout << "rk(B) = " << rk(B) << endl;
	cout << "det(A) = " << det(A) << endl;
	cout << endl;

	// operations
	C = A * B;
	C = 5 * C;
	D = B + C;
	E = A ^ -1; // inverse matrix

	// other functions
	F = steppedView(E);
	G = betterSteppedView(E);
	H = transpose(E);

	//////////////////////////////

	// printing matrix E
	cout << "Matrix E:\n" << endl;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
			cout << E(i, j) << "	";
		cout << endl;
	}
	cout << endl;

	//////////////////////////////

	for (int i = 0; i < 40; ++i)
		cout << "=";
	cout << "\n" << endl;

	//////////////////////////////
	/// VECTORS
	//////////////////////////////

	// constructing vector of a size 3
	vec k(3), l(3), m(3), p(3);

	// initializing vectors
	k(0) = 3.2; k(1) = 6.3; k(2) = 8.1;
	l(0) = 9.3; l(1) = 6.9; l(2) = 3.6;
	m(0) = 3.1; m(1) = 2.3; m(2) = 1.2;
	p = l + m;

	// functions
	cout << "Are l and p collinear? Answer = " << areCollinear(l, p) << endl;
	cout << "Are k, l and m coplanar? Answer =  " << areCoplanar(k, l, m) << endl;
	cout << "length(k) = " << length(k) << endl;
	cout << "cos between k and l = " << cosBetween(k, l) << endl;
	cout << "dot product of k and l = " << dotProd(k, l) << endl;
	cout << "cross product of k and l = " << crossProd(k, l) << endl;
	cout << endl;

	// turning vectors to matrices
	vec t(4);
	mat Q(1, 4);
	mat R(4, 1);
	Q = asRow(t);
	R = asCol(t);

	// converting matrices (of proper size) to vectors
	t = asVec(Q);
	t = asVec(R);

	//////////////////////////////

	for (int i = 0; i < 40; ++i)
		cout << "=";
	cout << "\n" << endl;

	//////////////////////////////
	/// COMPLEX NUMBERS
	//////////////////////////////

	// constructors
	complex z(3, 8);
	complex w, v, q;

	// ways of setting the complex number
	w = 4;
	v(5.5, 12.8);

	// real & imaginary parts
	cout << "v = " << v << endl;
	cout << "Re(v) = " << Re(v) << ", Im(v) = " << Im(v) << endl;
	cout << endl;

	// all operations
	q = z * w + v;
	cout << "q = " << q << endl;

	// module & argument
	cout << "mod(q) = " << mod(q) << endl;
	cout << "arg(q) = " << arg(q) << endl;

	// conjugated complex number
	cout << "conj(q) = " << conj(q) << endl;
	cout << endl;

	// complex roots of the power 3
	int n = 3;
	complex* roots = complRoots(q, n);

	// output
	cout << "Complex roots:\n" << endl;
	for (int i = 0; i < n; ++i)
		cout << roots[i] << endl;
	cout << endl;

	delete[] roots;

	//////////////////////////////

	for (int i = 0; i < 40; ++i)
		cout << "=";
	cout << "\n" << endl;
	
	//////////////////////////////
	/// INTERPOLATION
	//////////////////////////////

	// interpolation using Newton's formula
	double x[5], y[5];
	x[0] = 4.0; y[0] = 12.0;
	x[1] = 5.0; y[1] = 9.0;
	x[2] = 6.0; y[2] = 8.0;
	////// 7.0 must be 9.0
	x[3] = 8.0; y[3] = 12.0;
	x[4] = 9.0; y[4] = 17.0;
	double s = 7.0;

	// inputs, outputs, inputs-outputs array size, point to calculate the value at
	cout << "Interpolated value = " << getInterpolatedValue(x, y, 5, s) << endl;

	return 0;
}
