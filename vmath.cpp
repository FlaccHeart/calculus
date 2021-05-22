#include "vmath.h"

namespace vm
{
	bool areEqual(const double f, const double s)
	{
		const double small = 1e-7;

		if (abs(f) < small && abs(s) < small)
			return abs(f - s) < small;

		if (abs(f - s) <= eps * fmax(abs(f), abs(s)))
			return true;

		return false;
	}

	void restoreInt(double& n)
	{
		if (areEqual(n, round(n)))
			n = round(n);

		if (n == -0) n = 0;
	}

	int sgn(const double x)
	{
		if (areEqual(x, 0.0)) return 0;
		if (x > 0) return 1;
		return -1;
	}

	unsigned factorial(const unsigned n)
	{
		if (n == 0) return 1;

		unsigned k = 2, r = 1;

		while (k <= n)
		{
			r *= k;
			++k;
		}

		return r;
	}

	double interpolate(double* x, double* y, unsigned n, double in)
	{
		double* u = new double[n];
		u[0] = y[0];
		double fx = u[0], p = in - x[0];

		for (unsigned i = 1; i < n; ++i)
		{
			double fn = 0.0;
			double d = 1.0;

			for (unsigned j = 0; j < i; ++j)
			{
				double c = 1.0;

				for (unsigned k = 0; k < j; ++k)
					c *= x[i] - x[k];

				fn += u[j] * c;
				d *= x[i] - x[j];
			}

			u[i] = (y[i] - fn) / d;
			fx += u[i] * p;
			p *= in - x[i];
		}

		delete[] u;
		return fx;
	}

	void invSets(unsigned*** mem, const unsigned n)
	{
		if (n < 2) return;

		invSets(mem, n - 1);

		const unsigned f = factorial(n);
		const unsigned p = factorial(n - 1);

		for (unsigned i = 1; i < n; ++i) // clone S(n - 1) block n - 1 times
		{
			for (unsigned j = 0; j < p; ++j) // inside a block
			{
				unsigned u = 0;

				while (mem[j][u][0] || mem[j][u][1])
				{
					mem[i * p + j][u][0] = mem[j][u][0];
					mem[i * p + j][u][1] = mem[j][u][1];
					++u;
				}

				for (unsigned k = u; k < u + i; ++k)
				{
					mem[i * p + j][k][0] = (n - 1) - i + k - u;
					mem[i * p + j][k][1] = (n - 1);
				}
			}
		}
	}

	void trspGroup(unsigned** group, const unsigned n)
	{
		// group itself must have a size [n!] * [n]

		const unsigned f = factorial(n);
		const unsigned s = n * (n - 1) / 2;

		// creating a memory block of a size [n!] * [n * (n - 1) / 2] * [2]
		unsigned*** trsp = new unsigned** [f];

		for (unsigned i = 0; i < f; ++i)
		{
			trsp[i] = new unsigned* [s];

			for (unsigned j = 0; j < s; ++j)
			{
				trsp[i][j] = new unsigned[2];
				trsp[i][j][0] = trsp[i][j][1] = 0;
			}
		}

		// getting inversions to build transpositions
		invSets(trsp, n);

		for (unsigned i = 0; i < f; ++i)
		{
			// initializing a group element with id transposition
			for (unsigned j = 0; j < n; ++j)
				group[i][j] = j;

			// applying inversions
			for (unsigned u = 0; u < s; ++u)
			{
				if (!trsp[i][u][0] && !trsp[i][u][1])
					break;

				const unsigned t = group[i][trsp[i][u][0]];
				group[i][trsp[i][u][0]] = group[i][trsp[i][u][1]];
				group[i][trsp[i][u][1]] = t;
			}
		}

		for (unsigned i = 0; i < f; ++i)
		{
			for (unsigned j = 0; j < s; ++j)
				delete[] trsp[i][j];
			delete[] trsp[i];
		}

		delete[] trsp;
	}

	////////////////////////////////////////////////////
	/// MATRIX CLASS
	////////////////////////////////////////////////////

	// CONSTRUCTORS & DESTRUCTOR

	mat::mat()
	{
		rows = cols = 0;
		el = nullptr;
	}

	mat::mat(const unsigned _rows, const unsigned _cols)
	{
		if (_rows == 0 || _cols == 0)
		{
			rows = cols = 0;
			el = nullptr;
		}
		else
		{
			rows = _rows; cols = _cols;

			el = (double**)malloc(sizeof(double*) * rows + sizeof(double) * rows * cols);
			el[0] = (double*)(el + rows);

			for (unsigned j = 0; j < cols; ++j)
				el[0][j] = 0.0;

			for (unsigned i = 1; i < rows; ++i)
			{
				el[i] = el[i - 1] + cols;

				for (unsigned j = 0; j < cols; ++j)
					el[i][j] = 0.0;
			}
		}
	}

	mat::mat(const mat& o)
	{
		rows = o.rows; cols = o.cols;

		if (o.el == nullptr)
			el = nullptr;
		else
		{
			el = (double**)malloc(sizeof(double*) * rows + sizeof(double) * rows * cols);
			el[0] = (double*)(el + rows);

			for (unsigned j = 0; j < cols; ++j)
				el[0][j] = o.el[0][j];

			for (unsigned i = 1; i < rows; ++i)
			{
				el[i] = el[i - 1] + cols;

				for (unsigned j = 0; j < cols; ++j)
					el[i][j] = o.el[i][j];
			}
		}
	}

	mat::~mat()
	{
		free(el);
	}

	// OPERATORS

	// picking an element

	double& mat::operator()(const unsigned i, const unsigned j) const
	{
		if (i >= rows || j >= cols)
		{
			std::cout << "Out of range" << std::endl;
			exit(1);
		}

		return el[i][j];
	}

	// arithmetic

	mat mat::operator+(const mat& m) const
	{
		if (el == nullptr || m.el == nullptr)
		{
			std::cout << "Object is not defined" << std::endl;
			exit(1);
		}

		if (rows != m.rows || cols != m.cols)
		{
			std::cout << "Objects can not be added" << std::endl;
			exit(1);
		}

		mat sum(rows, cols);

		for (unsigned i = 0; i < rows; ++i)
			for (unsigned j = 0; j < cols; ++j)
				sum.el[i][j] = el[i][j] + m.el[i][j];

		return sum;
	}

	mat mat::operator-(const mat& m) const
	{
		if (el == nullptr || m.el == nullptr)
		{
			std::cout << "Object is not defined" << std::endl;
			exit(1);
		}

		if (rows != m.rows || cols != m.cols)
		{
			std::cout << "Objects can not be subtracted" << std::endl;
			exit(1);
		}

		mat dif(rows, cols);

		for (unsigned i = 0; i < rows; ++i)
			for (unsigned j = 0; j < cols; ++j)
				dif.el[i][j] = el[i][j] - m.el[i][j];

		return dif;
	}

	mat mat::operator-() const
	{
		return -1.0 * *this;
	}

	mat mat::operator*(const mat& m) const
	{
		if (el == nullptr || m.el == nullptr)
		{
			std::cout << "Object is not defined" << std::endl;
			exit(1);
		}

		if (cols != m.rows)
		{
			std::cout << "Matrices can not be multiplied" << std::endl;
			exit(1);
		}

		mat prod(rows, m.cols);

		for (unsigned i = 0; i < rows; ++i)
			for (unsigned j = 0; j < m.cols; ++j)
				for (unsigned k = 0; k < cols; ++k)
				{
					prod.el[i][j] += el[i][k] * m.el[k][j];
					restoreInt(prod.el[i][j]);
				}

		return prod;
	}

	// operations with constants

	mat operator*(const double lambda, const mat& m)
	{
		if (m.el == nullptr)
		{
			std::cout << "Object is not defined" << std::endl;
			exit(1);
		}

		mat res(m.rows, m.cols);

		for (unsigned i = 0; i < m.rows; ++i)
			for (unsigned j = 0; j < m.cols; ++j)
			{
				res.el[i][j] = lambda * m.el[i][j];
				restoreInt(res.el[i][j]);
			}

		return res;
	}

	mat operator*(const mat& m, const double lambda)
	{
		return operator*(lambda, m);
	}

	mat operator/(const mat& m, const double lambda)
	{
		if (m.el == nullptr)
		{
			std::cout << "Object is not defined" << std::endl;
			exit(1);
		}

		mat res(m.rows, m.cols);

		for (unsigned i = 0; i < m.rows; ++i)
			for (unsigned j = 0; j < m.cols; ++j)
			{
				res.el[i][j] = m.el[i][j] / lambda;
				restoreInt(res.el[i][j]);
			}

		return res;
	}

	mat mat::operator^(const int pow) const
	{
		if (el == nullptr)
		{
			std::cout << "Object is not defined" << std::endl;
			exit(1);
		}

		if (rows != cols)
		{
			std::cout << "Matrix can not be exponentiated" << std::endl;
			exit(1);
		}

		if (pow < -1 || pow == 0)
		{
			std::cout << "Powers can only be positive integers or -1 (for inverse)" << std::endl;
			exit(1);
		}

		mat res(rows, cols);

		if (pow > 0)
		{
			res = *this;
			for (int i = 1; i < pow; ++i)
				res *= *this;
		}
		else
		{
			if (rk(*this) != rows)
			{
				std::cout << "The inverse matrix does not exist" << std::endl;
				exit(1);
			}

			mat attr1(rows, 2 * cols), attr2(rows, 2 * cols);

			for (unsigned i = 0; i < rows; ++i)
			{
				for (unsigned j = 0; j < cols; ++j)
					attr1.el[i][j] = el[i][j];

				for (unsigned j = cols; j < 2 * cols; ++j)
					if (i == (j - cols))
						attr1.el[i][j] = 1.0;
			}

			attr2 = betterSteppedView(attr1);

			for (unsigned i = 0; i < rows; ++i)
				for (unsigned j = 0; j < cols; ++j)
					res(i, j) = attr2(i, j + cols);
		}
		return res;
	}

	// equating

	void mat::operator=(const mat& m)
	{
		if (el == nullptr) setSize(m.rows, m.cols);
		else if (rows != m.rows || cols != m.cols)
		{
			std::cout << "Objects can not be equated" << std::endl;
			exit(1);
		}

		for (unsigned i = 0; i < rows; ++i)
			for (unsigned j = 0; j < cols; ++j)
				el[i][j] = m.el[i][j];
	}

	void mat::operator+=(const mat& m)
	{
		*this = *this + m;
	}

	void mat::operator-=(const mat& m)
	{
		*this = *this - m;
	}

	void mat::operator*=(const mat& m)
	{
		*this = *this * m;
	}

	void mat::operator*=(const double lambda)
	{
		*this = *this * lambda;
	}

	void mat::operator/=(const double lambda)
	{
		*this = *this / lambda;
	}

	// comparison

	bool operator==(const mat& l, const mat& r)
	{
		if (l.rows != r.rows || l.cols != r.cols)
			return false;

		for (unsigned i = 0; i < l.rows; ++i)
			for (unsigned j = 0; j < l.cols; ++j)
				if (!areEqual(l.el[i][j], r.el[i][j]))
					return false;

		return true;
	}

	bool operator!=(const mat& l, const mat& r)
	{
		return !operator==(l, r);
	}

	// METHODS

	void mat::setSize(const unsigned m, const unsigned n)
	{
		mat t(m, n);

		for (unsigned i = 0; i < std::min(rows, m); ++i)
			for (unsigned j = 0; j < std::min(cols, n); ++j)
				t.el[i][j] = el[i][j];

		rows = m; cols = n;

		el = (double**)realloc(el, sizeof(double*) * rows + sizeof(double) * rows * cols);
		el[0] = (double*)(el + rows);

		for (unsigned j = 0; j < cols; ++j)
			el[0][j] = t.el[0][j];

		for (unsigned i = 1; i < rows; ++i)
		{
			el[i] = el[i - 1] + cols;

			for (unsigned j = 0; j < cols; ++j)
				el[i][j] = t.el[i][j];
		}
	}

	mat transpose(const mat& m)
	{
		if (m.el == nullptr)
		{
			std::cout << "Object is not defined" << std::endl;
			exit(1);
		}

		mat res(m.cols, m.rows);

		for (unsigned i = 0; i < m.cols; ++i)
			for (unsigned j = 0; j < m.rows; ++j)
				res(i, j) = m.el[j][i];

		return res;
	}

	mat steppedView(const mat& m)
	{
		if (m.el == nullptr)
		{
			std::cout << "Object is not defined" << std::endl;
			exit(1);
		}

		mat stepped(m.rows, m.cols);

		for (unsigned k = 0; k < m.rows - 1; ++k)
		{
			unsigned x = m.rows, y = m.cols;

			for (unsigned j = 0; j < m.cols; ++j)
			{
				for (unsigned i = k; i < m.rows; ++i)
					if (!areEqual(m.el[i][j], 0.0))
					{
						x = i; y = j;
						break;
					}

				if (x != m.rows) break;
			}

			stepped = m;

			if (x == m.rows) return stepped; // all zeros

			if (x != k)
				for (unsigned j = 0; j < m.cols; ++j)
					stepped(k, j) += m.el[x][j];

			for (unsigned i = k + 1; i < m.rows; ++i)
			{
				for (unsigned j = m.cols - 1; j > y; --j)
					stepped(i, j) -= m.el[i][y] / m.el[k][y] * m.el[k][j];
				stepped(i, y) -= m.el[i][y] / m.el[k][y] * m.el[k][y];
			}
		}

		for (unsigned i = 0; i < m.rows; ++i)
			for (unsigned j = 0; j < m.cols; ++j)
				restoreInt(stepped(i, j));

		return stepped;
	}

	mat betterSteppedView(const mat& m)
	{
		mat stepped = steppedView(m);
		unsigned rank = rk(stepped);

		for (unsigned i = 0; i < rank; ++i)
		{
			double lambda = 0.0;

			for (unsigned j = 0; j < stepped.cols; ++j)
			{
				if (!areEqual(stepped(i, j), 0.0) && areEqual(lambda, 0.0))
					lambda = stepped(i, j);

				if (!areEqual(lambda, 0.0))
					stepped(i, j) /= lambda;
			}
		}

		for (unsigned i = 0; i < rank; ++i)
		{
			unsigned y = stepped.cols;

			for (unsigned j = 0; j < stepped.cols; ++j)
				if (!areEqual(stepped(i, j), 0.0))
				{
					y = j;
					break;
				}

			for (unsigned p = 0; p < i; ++p)
			{
				double lambda = stepped(i, y) * stepped(p, y);

				for (unsigned q = y; q < stepped.cols; ++q)
					stepped(p, q) -= lambda * stepped(i, q);
			}
		}

		for (unsigned i = 0; i < stepped.rows; ++i)
			for (unsigned j = 0; j < stepped.cols; ++j)
				restoreInt(stepped(i, j));

		return stepped;
	}

	unsigned rk(const mat& m)
	{
		mat stepped = steppedView(m);
		unsigned q = 0;

		for (unsigned i = 0; i < stepped.rows; ++i)
		{
			bool found = false;

			for (unsigned j = 0; j < stepped.cols; ++j)
				if (!areEqual(stepped(i, j), 0.0))
				{
					found = true;
					break;
				}

			if (!found) return q;
			++q;
		}

		return q;
	}

	double det(const mat& m)
	{
		if (m.el == nullptr)
		{
			std::cout << "Object is not defined" << std::endl;
			exit(1);
		}

		if (m.rows != m.cols)
		{
			std::cout << "Determinant can only be calculated for square matrices" << std::endl;
			exit(1);
		}

		const unsigned n = m.rows;
		const unsigned f = factorial(n);
		unsigned** group = new unsigned* [f];

		for (unsigned i = 0; i < f; ++i)
			group[i] = new unsigned[n];

		trspGroup(group, n);

		double determinant = 0.0;

		for (unsigned i = 0; i < f; ++i)
		{
			int sign;

			if (i % 2 == 0)
			{
				if (i % 4 == 0) sign = 1;
				else sign = -1;
			}
			else
			{
				if ((i - 1) % 4 == 0) sign = -1;
				else sign = 1;
			}

			double term = 1.0;

			for (unsigned j = 0; j < n; ++j)
				term *= m.el[j][group[i][j]];

			term *= sign;
			determinant += term;
		}

		for (unsigned i = 0; i < f; ++i)
			delete[] group[i];
		delete[] group;

		restoreInt(determinant);
		return determinant;
	}

	////////////////////////////////////////////////////
	/// VECTOR CLASS
	////////////////////////////////////////////////////

	// CONSTRUCTORS & DESTRUCTOR

	vec::vec()
	{
	}

	vec::vec(unsigned n)
	{
		m.setSize(n, 1);
	}

	vec::vec(const vec& o)
	{
		m.setSize(o.m.rows, 1);
		m = o.m;
	}

	vec::~vec()
	{
	}

	// OPERATORS

	// picking an element

	double& vec::operator[](const unsigned k) const
	{
		return m.el[k][0];
	}

	// arithmetic

	vec vec::operator+(const vec& v) const
	{
		vec sum(v.m.rows);
		sum.m = this->m + v.m;
		return sum;
	}

	vec vec::operator-(const vec& v) const
	{
		vec dif(v.m.rows);
		dif.m = this->m - v.m;
		return dif;
	}

	vec vec::operator-() const
	{
		vec res(m.rows);
		res = -1 * m;
		return res;
	}

	// operations with constants

	vec operator*(const double lambda, const vec& v)
	{
		vec res(v.m.rows);
		res.m = lambda * v.m;
		return res;
	}

	vec operator*(const vec& v, const double lambda)
	{
		return operator*(lambda, v);
	}

	vec operator/(const vec& v, const double lambda)
	{
		vec res(v.m.rows);
		res.m = v.m / lambda;
		return res;
	}

	// equating

	void vec::operator=(const vec& v)
	{
		m = v.m;
	}

	void vec::operator=(const mat& mt)
	{
		if (mt.cols == 1)
			m = mt;
		else if (mt.rows == 1)
			m = transpose(mt);
		else
		{
			std::cout << "This matrix can not be converted to vector" << std::endl;
			exit(1);
		}
	}

	void vec::operator+=(const vec& v)
	{
		*this = *this + v;
	}

	void vec::operator-=(const vec& v)
	{
		*this = *this - v;
	}

	void vec::operator*=(const double lambda)
	{
		*this = *this * lambda;
	}

	void vec::operator/=(const double lambda)
	{
		*this = *this / lambda;
	}

	// comparison

	bool operator==(const vec& l, const vec& r)
	{
		return l.m == r.m;
	}

	bool operator!=(const vec& l, const vec& r)
	{
		return l.m != r.m;
	}

	// ostream

	std::ostream& operator<<(std::ostream& out, const vec& v)
	{
		out << "{ ";

		for (unsigned i = 0; i < v.m.rows - 1; ++i)
			out << v.m.el[0][i] << ", ";
		out << v.m.el[0][v.m.rows - 1] << " }";

		return out;
	}

	// METHODS

	mat asRow(const vec& v)
	{
		return transpose(v.m);
	}

	mat asCol(const vec& v)
	{
		return v.m;
	}

	void vec::setSize(const unsigned n)
	{
		m.setSize(n, 1);
	}

	bool areCollinear(const vec& a, const vec& b)
	{
		if (a.m.rows != b.m.rows)
		{
			std::cout << "Ñollinearity is not defined for vectors of different dimensions" << std::endl;
			exit(1);
		}

		bool nonZeroA = false, nonZeroB = false;

		for (unsigned i = 0; i < a.m.rows; ++i)
		{
			if (!areEqual(a.m.el[i][0], 0.0)) nonZeroA = true;
			if (!areEqual(b.m.el[i][0], 0.0)) nonZeroB = true;
		}

		if (!nonZeroA || !nonZeroB)
			return true; // both are null vectors

		double k = 0.0;

		for (unsigned i = 0; i < a.m.rows; ++i)
		{
			if (!areEqual(b.m.el[i][0], 0.0))
			{
				if (areEqual(a.m.el[i][0], 0.0)) return false;
				else
				{
					if (!areEqual(k, 0.0))
					{
						if (!areEqual(a.m.el[i][0] / b.m.el[i][0], k))
							return false;
					}
					else
						k = a.m.el[i][0] / b.m.el[i][0];
				}
			}
			else if (!areEqual(a.m.el[i][0], 0.0)) return false;
		}

		return true;
	}

	bool areCoplanar(const vec& a, const vec& b, const vec& c)
	{
		if (a.m.rows != b.m.rows || a.m.rows != c.m.rows || a.m.rows != 3)
		{
			std::cout << "bool areCoplanar(...) is only defined for 3 dimensional vectors" << std::endl;
			exit(1);
		}

		mat D(3, 3);

		D(0, 0) = a.m.el[0][0]; D(0, 1) = a.m.el[1][0]; D(0, 2) = a.m.el[2][0];
		D(1, 0) = b.m.el[0][0]; D(1, 1) = b.m.el[1][0]; D(1, 2) = b.m.el[2][0];
		D(2, 0) = c.m.el[0][0]; D(2, 1) = c.m.el[1][0]; D(2, 2) = c.m.el[2][0];

		return !det(D);
	}

	double length(const vec& v)
	{
		double m = 0.0;

		for (unsigned i = 0; i < v.m.rows; ++i)
			m += pow(v.m.el[i][0], 2);

		return pow(m, 0.5);
	}

	double dotProd(const vec& a, const vec& b)
	{
		return (a.m * transpose(b.m)).el[0][0];
	}

	double cosBetween(const vec& a, const vec& b)
	{
		return dotProd(a, b) / length(a) / length(b);
	}

	vec crossProd(const vec& a, const vec& b)
	{
		if (a.m.rows != 3)
		{
			std::cout << "Cross product is only defined for 3 dimensional vectors" << std::endl;
			exit(1);
		}

		vec res(3);

		mat x(2, 2), y(2, 2), z(2, 2);

		x(0, 0) = a.m.el[1][0]; x(0, 1) = a.m.el[2][0];
		x(1, 0) = b.m.el[1][0]; x(1, 1) = b.m.el[2][0];

		y(0, 0) = a.m.el[2][0]; y(0, 1) = a.m.el[0][0];
		y(1, 0) = b.m.el[2][0]; y(1, 1) = b.m.el[0][0];

		z(0, 0) = a.m.el[0][0]; z(0, 1) = a.m.el[1][0];
		z(1, 0) = b.m.el[0][0]; z(1, 1) = b.m.el[1][0];

		res[0] = det(x); res[1] = det(y); res[2] = det(z);

		return res;
	}

	////////////////////////////////////////////////////
	/// COMPLEX CLASS
	////////////////////////////////////////////////////

	// CONSTRUCTORS & DESTRUCTOR

	complex::complex()
	{
		m.setSize(2, 2);
		m(0, 0) = 0.0; m(0, 1) = 0.0;
		m(1, 0) = 0.0; m(1, 1) = 0.0;
	}

	complex::complex(const double re, const double im)
	{
		m.setSize(2, 2);
		m(0, 0) = re; m(0, 1) = im;
		m(1, 0) = -im; m(1, 1) = re;
	}

	complex::complex(const complex& o)
	{
		m.setSize(o.m.rows, o.m.cols);
		m = o.m;
	}

	complex::~complex()
	{
	}

	// OPERATORS

	// arithmetic

	void complex::operator()(const double re, const double im) const
	{
		m(0, 0) = re; m(0, 1) = im;
		m(1, 0) = -im; m(1, 1) = re;
	}

	complex complex::operator+(const complex& c) const
	{
		complex sum;
		sum.m = m + c.m;
		return sum;
	}

	complex complex::operator-(const complex& c) const
	{
		complex dif;
		dif.m = m - c.m;
		return dif;
	}

	complex complex::operator*(const complex& c) const
	{
		complex prod;
		prod.m = m * c.m;
		return prod;
	}

	complex complex::operator+(const double re) const
	{
		return *this + complex(re, 0.0);
	}

	complex complex::operator-(const double re) const
	{
		return *this - complex(re, 0.0);
	}

	complex complex::operator*(const double re) const
	{
		return *this * complex(re, 0.0);
	}

	// operators with constant

	complex complex::operator^(const int pow) const
	{
		complex res;
		res.m = m ^ pow;
		return res;
	}

	// equating

	void complex::operator=(const complex& c)
	{
		m = c.m;
	}

	void complex::operator+=(const complex& c)
	{
		*this = *this + c;
	}

	void complex::operator-=(const complex& c)
	{
		*this = *this - c;
	}

	void complex::operator*=(const complex& c)
	{
		*this = *this * c;
	}

	void complex::operator=(const double re)
	{
		m(0, 0) = re; m(0, 1) = 0.0;
		m(1, 0) = 0.0; m(1, 1) = re;
	}

	void complex::operator+=(const double re)
	{
		*this = *this + complex(re, 0.0);
	}

	void complex::operator-=(const double re)
	{
		*this = *this - complex(re, 0.0);
	}

	void complex::operator*=(const double re)
	{
		*this = *this * complex(re, 0.0);
	}

	// comparison

	bool operator==(const complex& l, const complex& r)
	{
		return l.m == r.m;
	}

	bool operator!=(const complex& l, const complex& r)
	{
		return l.m != r.m;
	}

	// ostream

	std::ostream& operator<<(std::ostream& out, const complex& c)
	{
		out << c.m.el[0][0] << " + " << c.m.el[0][1] << "i";
		return out;
	}

	// METHODS

	void complex::setRe(const double re)
	{
		m(0, 0) = re;
		m(1, 1) = re;
	}

	void complex::setIm(const double im)
	{
		m(0, 1) = im;
		m(1, 0) = -im;
	}

	double Re(const complex& c)
	{
		return c.m.el[0][0];
	}

	double Im(const complex& c)
	{
		return c.m.el[0][1];
	}

	double mod(const complex& c)
	{
		return pow(pow(Re(c), 2) + pow(Im(c), 2), 0.5);
	}

	double arg(const complex& c)
	{
		return acos(Re(c) / mod(c));
	}

	complex conj(const complex& c)
	{
		complex conjugated;
		conjugated.m = transpose(c.m);
		return conjugated;
	}

	complex* complRoots(const complex& c, const int p)
	{
		if (p <= 0)
		{
			std::cout << "You can only take the complex roots of the positive powers" << std::endl;
			exit(1);
		}

		complex* roots = new complex[p];
		double newMod = pow(mod(c), 1 / (double)p);
		double newArg;

		for (int i = 0; i < p; ++i)
		{
			newArg = (arg(c) + PI * i) / (double)p;
			roots[i].m.el[0][0] = newMod * cos(newArg);
			roots[i].m.el[0][1] = newMod * sin(newArg);
			roots[i].m.el[1][0] = -newMod * sin(newArg);
			roots[i].m.el[1][1] = newMod * cos(newArg);

			for (unsigned j = 0; j < roots[i].m.rows; ++j)
				for (unsigned k = 0; k < roots[i].m.cols; ++k)
					restoreInt(roots[i].m.el[j][k]);
		}

		return roots;
	}
}
