#include <iostream>
#include "vmath.h"

using namespace std;

int main()
{
	//////////////////////////////
	/// MATRICES
	//////////////////////////////

	// matrices constructor examples
	vm::mat A(3, 3), B(3, 2), C(3, 2), D(3, 3), E(3, 3), F(3, 3), G(3, 3), H(3, 3);

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
	C += B;
	D = A ^ -1; // inverse matrix

	// other functions
	E = steppedView(D);
	F = betterSteppedView(D);
	G = transpose(D);
	D.setSize(3, 4); // fills new column with zeroes

	//////////////////////////////

	// printing matrix E
	cout << "Matrix D:\n" << endl;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 4; ++j)
			cout << D(i, j) << "	";
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
	vm::vec k(3), l(3), m(3), p(3);

	// initializing vectors
	k[0] = 3.2; k[1] = 6.3; k[2] = 8.1;
	l[0] = 9.3; l[1] = 6.9; l[2] = 3.6;
	m[0] = 3.1; m[1] = 2.3; m[2] = 1.2;
	p = l + m;

	// functions
	cout << "Are l and p collinear? Answer = " << areCollinear(l, p) << endl;
	cout << "Are k, l and m coplanar? Answer = " << areCoplanar(k, l, m) << endl;
	cout << "length(k) = " << length(k) << endl;
	cout << "cos between k and l = " << cosBetween(k, l) << endl;
	cout << "dot product of k and l = " << dotProd(k, l) << endl;
	cout << "cross product of k and l = " << crossProd(k, l) << endl;
	cout << endl;

	// turning vectors to matrices
	vm::vec t(4);
	t[0] = 4.65;
	t[1] = 2.92;
	t[2] = 8.21;
	t[3] = 5.27;
	vm::mat Q(1, 4);
	vm::mat R(4, 1);
	Q = asRow(t);
	R = asCol(t);

	// converting 1*n or n*1 matrices to vectors
	t = Q;
	t = R;

	// changing size
	cout << t << endl;
	t.setSize(3);
	cout << t << endl;
	t.setSize(4);
	cout << t << endl;
	cout << endl;

	//////////////////////////////

	for (int i = 0; i < 40; ++i)
		cout << "=";
	cout << "\n" << endl;

	//////////////////////////////
	/// COMPLEX NUMBERS
	//////////////////////////////

	// constructors
	vm::complex z(3, 8); // Re(w) == 3.0, Im(w) == 8.0
	vm::complex w, v, q;

	// ways of setting the complex number
	w = 4.1; // Re(w) == 4.1, Im(w) == 0.0
	w.setRe(3.9); // changing real part (only)
	w.setIm(7.3); // changing imaginary part (only)
	v(5.5, 12.8); // setting Re and Im parts simultaneously

	// real & imaginary parts
	cout << "v = " << v << endl;
	cout << "Re(v) = " << Re(v) << ", Im(v) = " << Im(v) << endl;
	cout << endl;

	// all operations
	q = z * w + v;
	cout << "q = " << q << endl; // "<<" is overloaded

	// module & argument
	cout << "mod(q) = " << mod(q) << endl;
	cout << "arg(q) = " << arg(q) << endl;

	// conjugated complex number
	cout << "conj(q) = " << conj(q) << endl;
	cout << endl;

	// complex roots of the power 3
	int n = 3;
	vm::complex* roots = complRoots(q, n);

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

	// inputs, outputs, in-out array size, point to calculate the value at
	cout << "Interpolated value = " << vm::getInterpolatedValue(x, y, 5, s) << endl;

	return 0;
}
