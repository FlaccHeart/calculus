#include "calculus.h"
#define eps 1e-15

// double comparison
bool comp_d(const double f, const double s)
{
	double small = 1e-7;
	if (fabs(f) < small && fabs(s) < small)
		return fabs(f - s) < small;

	if (fabs(f - s) <= eps * fmax(fabs(f), fabs(s)))
		return true;
	return false;
}

////////////////////////////////////////////////////
/// MATRIX CLASS
////////////////////////////////////////////////////

// default constructor

mat::mat(const unsigned _rows, const unsigned _cols)
{
	rows = _rows; cols = _cols;

	el = new double* [rows];
	el[0] = new double[rows * cols];

	for (unsigned j = 0; j < cols; ++j)
		el[0][j] = 0.0;

	for (unsigned i = 1; i < rows; ++i)
	{
		el[i] = el[0] + i * cols;

		for (unsigned j = 0; j < cols; ++j)
			el[i][j] = 0.0;
	}
}

// copy constructor

mat::mat(const mat& o)
{
	rows = o.rows; cols = o.cols;
	el = new double* [rows];
	el[0] = new double[rows * cols];

	for (unsigned j = 0; j < cols; ++j)
		el[0][j] = o.el[0][j];

	for (unsigned i = 1; i < rows; ++i)
	{
		el[i] = el[0] + i * cols;

		for (unsigned j = 0; j < cols; ++j)
			el[i][j] = o.el[i][j];
	}
}

// destructor

mat::~mat()
{
	delete[] el[0];
	delete[] el;
}

// operators overloading

double& mat::operator()(const unsigned i, const unsigned j)
{
	return el[i][j];
}

mat mat::operator+(const mat& m) const
{
	if (rows != m.rows || cols != m.cols)
	{
		std::cout << "Matrices cannot be added" << std::endl;
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
	if (rows != m.rows || cols != m.cols)
	{
		std::cout << "Matrices cannot be added" << std::endl;
		exit(1);
	}

	mat dif(rows, cols);

	for (unsigned i = 0; i < rows; ++i)
		for (unsigned j = 0; j < cols; ++j)
			dif.el[i][j] = el[i][j] - m.el[i][j];

	return dif;
}

mat mat::operator*(const mat& m) const
{
	if (cols != m.rows)
	{
		std::cout << "Matrices cannot be multiplied" << std::endl;
		exit(1);
	}

	mat prod(rows, m.cols);

	for (unsigned i = 0; i < rows; ++i)
		for (unsigned j = 0; j < m.cols; ++j)
			for (unsigned k = 0; k < cols; ++k)
				prod.el[i][j] += el[i][k] * m.el[k][j];

	return roundByEps(prod);
}

mat operator*(const double lambda, const mat& m)
{
	mat res(m.rows, m.cols);

	for (unsigned i = 0; i < m.rows; ++i)
		for (unsigned j = 0; j < m.cols; ++j)
			res.el[i][j] = lambda * m.el[i][j];

	return roundByEps(res);
}

mat operator*(const mat& m, const double lambda)
{
	return operator*(lambda, m);
}

mat mat::operator^(const int pow) const
{
	if (rows != cols)
	{
		std::cout << "Matrix cannot be exponentiated" << std::endl;
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
			res = res * *this;
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
				if (i == j - cols)
					attr1.el[i][j] = 1.0;
		}
		attr2 = betterSteppedView(attr1);

		for (unsigned i = 0; i < rows; ++i)
			for (unsigned j = 0; j < cols; ++j)
				res(i, j) = attr2(i, j + cols);
	}
	return roundByEps(res);
}

void mat::operator=(const mat& m)
{
	if (rows != m.rows || cols != m.cols)
	{
		std::cout << "Matrices cannot be equated" << std::endl;
		exit(1);
	}

	for (unsigned i = 0; i < rows; ++i)
		for (unsigned j = 0; j < cols; ++j)
			el[i][j] = m.el[i][j];
}

// friends

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

mat steppedView(mat m)
{
	for (unsigned k = 0; k < m.rows - 1; ++k)
	{
		unsigned x = m.rows, y = m.cols;
		for (unsigned j = 0; j < m.cols; ++j)
		{
			for (unsigned i = k; i < m.rows; ++i)
				if (!comp_d(m.el[i][j], 0.0))
				{
					x = i; y = j;
					break;
				}
			if (x != m.rows) break;
		}
		if (x == m.rows) return m; // all zeros

		if (x != k)
			for (unsigned j = 0; j < m.cols; ++j)
				m.el[k][j] += m.el[x][j];

		for (unsigned i = k + 1; i < m.rows; ++i)
		{
			for (unsigned j = m.cols - 1; j > y; --j)
				m.el[i][j] -= m.el[i][y] / m.el[k][y] * m.el[k][j];
			m.el[i][y] -= m.el[i][y] / m.el[k][y] * m.el[k][y];
		}
	}
	return roundByEps(m);
}

mat betterSteppedView(const mat& m)
{
	mat stepped(m.rows, m.cols);
	stepped = steppedView(m);
	int nonZeroQ = stepped.countNonZero(stepped);

	for (int i = 0; i < nonZeroQ; ++i)
	{
		double lambda = 0.0;
		for (unsigned j = 0; j < stepped.cols; ++j)
		{
			if (!comp_d(stepped(i, j), 0.0) && comp_d(lambda, 0.0))
				lambda = stepped(i, j);
			if (lambda != 0)
				stepped(i, j) /= lambda;
		}
	}

	for (int i = 0; i < nonZeroQ; ++i)
	{
		unsigned y = stepped.cols;
		for (unsigned j = 0; j < stepped.cols; ++j)
			if (!comp_d(stepped(i, j), 0.0))
			{
				y = j;
				break;
			}
		for (int p = 0; p < i; ++p)
		{
			double lambda = stepped(i, y) * stepped(p, y);
			for (unsigned q = y; q < stepped.cols; ++q)
				stepped(p, q) -= lambda * stepped(i, q);
		}
	}

	return roundByEps(stepped);
}

mat transpose(const mat& m)
{
	if (m.rows != m.cols)
	{
		std::cout << "You can only transpose square matrices" << std::endl;
		exit(1);
	}
	mat res(m.rows, m.cols);

	for (unsigned i = 0; i < m.rows; ++i)
		for (unsigned j = 0; j < m.cols; ++j)
			res(i, j) = m.el[j][i];

	return res;
}

int rk(mat m)
{
	return m.countNonZero(steppedView(m));
}

double det(mat m)
{
	if (m.rows != m.cols)
	{
		std::cout << "Determinant can only be calculated for square matrices" << std::endl;
		exit(1);
	}

	const unsigned n = m.rows;
	double determinant = 0.0;
	unsigned f = factorial(n);
	unsigned** p = m.transpositions();

	// calculating the determinant
	for (unsigned i = 0; i < f; ++i)
	{
		int sign;
		if (i % 2 == 0)
		{
			if (i % 4 == 0)
				sign = 1;
			else
				sign = -1;
		}
		else
		{
			if ((i - 1) % 4 == 0)
				sign = -1;
			else
				sign = 1;
		}

		double term = 1;
		for (unsigned j = 0; j < n; ++j)
			term *= m(j, p[i][j]);
		term *= sign;
		determinant += term;
	}

	delete[] p[0];
	delete[] p;

	return determinant;
}

int mat::countNonZero(mat m)
{
	int q = 0;
	for (unsigned i = 0; i < m.rows; ++i)
	{
		bool found = false;
		for (unsigned j = 0; j < m.cols; ++j)
			if (!comp_d(m(i, j), 0.0))
			{
				found = true;
				break;
			}
		if (!found) return q;
		++q;
	}
	return q;
}

void mat::swap(unsigned& f, unsigned& s)
{
	unsigned t = f;
	f = s;
	s = t;
}

unsigned** mat::transpositions()
{
	const unsigned n = rows;
	unsigned f = factorial(n);

	// inversions matrix (3rd dimension for transposition itself)
	unsigned*** invSets = new unsigned** [f];
	for (unsigned i = 0; i < f; ++i)
	{
		invSets[i] = new unsigned* [n * (n - 1) / 2];
		for (unsigned j = 0; j < n * (n - 1) / 2; ++j)
		{
			invSets[i][j] = new unsigned[2];
			invSets[i][j][0] = 0;
			invSets[i][j][1] = 0;
		}
	}

	// creating the sets of inversions

	invSets[1][0][0] = n - 2;
	invSets[1][0][1] = n - 1;

	for (unsigned i = 3; i <= n; ++i)
	{
		for (unsigned p = factorial(i - 1); p < factorial(i); ++p)
		{
			for (unsigned q = 0; q < i; ++q)
			{
				invSets[p][q][0] = n - i;
				invSets[p][q][1] = n - i + q + 1;
			}
		}

		for (unsigned j = 1; j < i; ++j)
		{
			for (unsigned k = 0; k < factorial(i - 1); ++k)
			{
				unsigned inv = n * (n - 1) / 2;

				for (unsigned r = 0; r < n * (n - 1) / 2; ++r)
					if (invSets[k][r][0] == invSets[k][r][1])
					{
						inv = r;
						break;
					}

				for (unsigned l = 0; l < inv; ++l)
				{
					invSets[factorial(i - 1) * j + k][j + l][0] = invSets[k][l][0];
					invSets[factorial(i - 1) * j + k][j + l][1] = invSets[k][l][1];
				}

				for (unsigned s = j + inv; s < n * (n - 1) / 2; ++s)
				{
					invSets[factorial(i - 1) * j + k][s][0] = 0;
					invSets[factorial(i - 1) * j + k][s][1] = 0;
				}
			}
		}
	}

	// creating transpositions matrix

	unsigned** trsp = new unsigned* [f];
	trsp[0] = new unsigned[n * f];

	for (unsigned j = 0; j < n; ++j)
		trsp[0][j] = j;

	for (unsigned i = 1; i < f; ++i)
	{
		trsp[i] = trsp[0] + i * n;
		for (unsigned j = 0; j < n; ++j)
			trsp[i][j] = j; // filling everything with "e" transpositions
	}

	// applying the inversions to transpositions

	for (unsigned i = 0; i < f; ++i)
		for (unsigned j = 0; j < n * (n - 1) / 2; ++j)
			swap(trsp[i][invSets[i][j][0]], trsp[i][invSets[i][j][1]]);

	for (unsigned i = 0; i < f; ++i)
	{
		for (unsigned j = 0; j < n * (n - 1) / 2; ++j)
			delete[] invSets[i][j];
		delete[] invSets[i];
	}
	delete[] invSets;

	return trsp;
}

mat roundByEps(const mat& m)
{
	for (unsigned i = 0; i < m.rows; ++i)
		for (unsigned j = 0; j < m.cols; ++j)
		{
			if (comp_d(m.el[i][j], round(m.el[i][j])))
				m.el[i][j] = round(m.el[i][j]);
			if (m.el[i][j] == -0)
				m.el[i][j] = 0;
		}
	return m;
}

////////////////////////////////////////////////////
/// COMPLEX CLASS
////////////////////////////////////////////////////

complex::complex()
{
	m = new mat(2, 2);
	(*m)(0, 0) = 0.0; (*m)(0, 1) = 0.0;
	(*m)(1, 0) = 0.0; (*m)(1, 1) = 0.0;
}

complex::complex(const double re, const double im)
{
	m = new mat(2, 2);
	(*m)(0, 0) = re; (*m)(0, 1) = im;
	(*m)(1, 0) = -im; (*m)(1, 1) = re;
}

complex::complex(const complex& o)
{
	m = new mat(2, 2);
	*m = *(o.m);
}

complex::~complex()
{
	delete m;
}

void complex::operator()(const double re, const double im)
{
	(*m)(0, 0) = re; (*m)(0, 1) = im;
	(*m)(1, 0) = -im; (*m)(1, 1) = re;
}

complex complex::operator+(const complex& c) const
{
	complex sum;
	*(sum.m) = *m + *(c.m);
	return sum;
}

complex complex::operator-(const complex& c) const
{
	complex dif;
	*(dif.m) = *m - *(c.m);
	return dif;
}

complex complex::operator*(const complex& c) const
{
	complex prod;
	*(prod.m) = *m * *(c.m);
	return prod;
}

complex complex::operator^(const int pow) const
{
	complex res;
	*(res.m) = *m ^ pow;
	return res;
}

void complex::operator=(const complex& c)
{
	*m = *(c.m);
}

void complex::operator=(const double re)
{
	(*m)(0, 0) = re; (*m)(0, 1) = 0.0;
	(*m)(1, 0) = 0.0; (*m)(1, 1) = re;
}

std::ostream& operator<<(std::ostream& out, const complex& c)
{
	out << c.m->operator()(0, 0) << " + " << c.m->operator()(0, 1) << "i";
	return out;
}

double& Re(const complex& c)
{
	return c.m->operator()(0, 0);
}

double& Im(const complex& c)
{
	return c.m->operator()(0, 1);
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
	*(conjugated.m) = transpose(*(c.m));
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
		roots[i].m->operator()(0, 0) = newMod * cos(newArg);
		roots[i].m->operator()(0, 1) = newMod * sin(newArg);
		roots[i].m->operator()(1, 0) = -newMod * sin(newArg);
		roots[i].m->operator()(1, 1) = newMod * cos(newArg);
		*(roots[i].m) = roundByEps(*(roots[i].m));
	}
	return roots;
}

double getInterpolatedValue(double* x, double* y, unsigned n, double in)
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
	return fx;
}
