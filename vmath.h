#pragma once
#include <iostream>
#include <cmath>

namespace vm
{
	const double PI = 3.141592653589793;
	const double eps = 1e-15;

	class mat;
	class vec;
	class complex;

	bool areEqual(const double, const double);
	void restoreInt(double&);
	int sgn(const double);
	unsigned factorial(const unsigned);
	double interpolate(double*, double*, unsigned, double);
	void invSets(unsigned***, const unsigned);
	void trspGroup(unsigned**, const unsigned);

	class mat
	{
	public:

		mat();
		mat(const unsigned, const unsigned);
		mat(const mat&);
		~mat();

		unsigned rows, cols;

		// picking an element
		double& operator()(const unsigned, const unsigned) const;

		// arithmetic
		mat operator+(const mat&) const;
		mat operator-(const mat&) const;
		mat operator-() const;
		mat operator*(const mat&) const;

		// operations with constants
		friend mat operator*(const double, const mat&);
		friend mat operator*(const mat&, const double);
		friend mat operator/(const mat&, const double);
		mat operator^(const int) const;

		// equating
		void operator=(const mat&);
		void operator+=(const mat&);
		void operator-=(const mat&);
		void operator*=(const mat&);
		void operator*=(const double);
		void operator/=(const double);

		// comparison
		friend bool operator==(const mat&, const mat&);
		friend bool operator!=(const mat&, const mat&);

		// methods
		void setSize(const unsigned, const unsigned);
		friend mat transpose(const mat&);
		friend mat steppedView(const mat&);
		friend mat betterSteppedView(const mat&);
		friend unsigned rk(const mat&);
		friend double det(const mat&);

		// elements
		double** el;
	};

	class vec
	{
	public:

		// constructors & destructor
		vec();
		vec(const unsigned);
		vec(const vec&);
		~vec();

		// picking an element
		double& operator[](const unsigned) const;

		// arithmetic
		vec operator+(const vec&) const;
		vec operator-(const vec&) const;
		vec operator-() const;

		// operations with constants
		friend vec operator*(const double, const vec&);
		friend vec operator*(const vec&, const double);
		friend vec operator/(const vec&, const double);

		// equating
		void operator=(const vec&);
		void operator=(const mat&);
		void operator+=(const vec&);
		void operator-=(const vec&);
		void operator*=(const double);
		void operator/=(const double);

		// comparison
		friend bool operator==(const vec&, const vec&);
		friend bool operator!=(const vec&, const vec&);

		// ostream
		friend std::ostream& operator<<(std::ostream&, const vec&);

		// methods
		friend mat asRow(const vec&);
		friend mat asCol(const vec&);
		void setSize(const unsigned);
		friend bool areCollinear(const vec&, const vec&);
		friend bool areCoplanar(const vec&, const vec&, const vec&);
		friend double length(const vec&);
		friend double dotProd(const vec&, const vec&);
		friend double cosBetween(const vec&, const vec&);
		friend vec crossProd(const vec&, const vec&);

	private:

		mat m;
	};

	class complex
	{
	public:

		// constructors & destructor
		complex();
		complex(const double, const double);
		complex(const complex&);
		~complex();

		// setting the number
		void operator()(const double, const double) const;

		// arithmetic
		complex operator+(const complex&) const;
		complex operator-(const complex&) const;
		complex operator*(const complex&) const;
		complex operator+(const double) const;
		complex operator-(const double) const;
		complex operator*(const double) const;

		// operators with constant
		complex operator^(const int) const;

		// equating
		void operator=(const complex&);
		void operator+=(const complex&);
		void operator-=(const complex&);
		void operator*=(const complex&);
		void operator=(const double);
		void operator+=(const double);
		void operator-=(const double);
		void operator*=(const double);

		// comparison
		friend bool operator==(const complex&, const complex&);
		friend bool operator!=(const complex&, const complex&);

		// ostream
		friend std::ostream& operator<<(std::ostream&, const complex&);

		// methods
		void setRe(const double);
		void setIm(const double);
		friend double Re(const complex&);
		friend double Im(const complex&);
		friend double mod(const complex&);
		friend double arg(const complex&);
		friend complex conj(const complex&);
		friend complex* complRoots(const complex&, const int);

	private:

		mat m;
	};
}
