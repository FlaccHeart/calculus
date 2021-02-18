#pragma once
#include <iostream>
#include <cmath>

namespace vm
{
	const double PI = 3.141592653589793;

	class mat;
	class vec;
	class complex;

	// functions
	bool areEqual(const double, const double);
	void roundByEps(double&);
	int sgn(double);
	unsigned factorial(const unsigned);
	double getInterpolatedValue(double*, double*, unsigned, double);

	class mat
	{
	public:

		// constructors & destructor
		mat();
		mat(const unsigned, const unsigned);
		mat(const mat&);
		~mat();

		unsigned rows, cols;
				
		// picking the element
		double& operator()(const unsigned, const unsigned);

		// arithmetic
		mat operator+(const mat&) const;
		mat operator-(const mat&) const;
		mat operator-() const;
		mat operator*(const mat&) const;

		// with constants
		friend mat operator*(const double, const mat&);
		friend mat operator*(const mat&, const double);
		friend mat operator/(const mat&, const double);
		mat operator^(const int) const;

		// equating
		void operator=(const mat&);
		void operator+=(const mat&);
		void operator-=(const mat&);
		void operator*=(const double);
		void operator/=(const double);

		// comparison
		friend bool operator==(const mat&, const mat&);
		friend bool operator!=(const mat&, const mat&);

		// methods
		void setSize(const unsigned, const unsigned);
		void swap(unsigned&, unsigned&);
		friend mat transpose(const mat&);
		friend mat steppedView(mat);
		friend mat betterSteppedView(const mat&);
		friend int rk(mat);
		friend double det(mat);

		// container
		double** el;

	private:

		unsigned** transpositions();
		int countNonZero(mat);
	};

	class vec
	{
	public:

		// constructors & destructor
		vec();
		vec(const unsigned);
		vec(const vec&);
		~vec();

		// picking the element
		double& operator[](const unsigned);

		// arithmetic
		vec operator+(const vec&) const;
		vec operator-(const vec&) const;
		vec operator-();

		// with constants
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
		friend double cosBetween(const vec&, const vec&);
		friend double dotProd(const vec&, const vec&);
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
		void operator()(const double, const double);

		// arithmetic
		complex operator+(const complex&) const;
		complex operator+(const double) const;
		complex operator-(const complex&) const;
		complex operator-(const double) const;
		complex operator*(const complex&) const;
		complex operator*(const double) const;

		// with constant
		complex operator^(const int) const;

		// equating
		void operator=(const complex&);
		void operator=(const double);
		void operator+=(const complex&);
		void operator+=(const double);
		void operator-=(const complex&);
		void operator-=(const double);
		void operator*=(const complex&);
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

} // namespace vm
