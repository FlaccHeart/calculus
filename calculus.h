#pragma once
#include <iostream>
#include <cmath>

const double PI = 3.141592653589793;
bool isEqual(const double, const double);
void roundByEps(double&);

class mat;
class vec;
class complex;

class mat
{
public:

	mat(const unsigned, const unsigned);
	mat(const mat&);
	~mat();

	unsigned rows, cols;

	// operators overloading
	double& operator()(const unsigned, const unsigned);
	mat operator+(const mat&) const;
	mat operator-(const mat&) const;
	mat operator*(const mat&) const;
	friend mat operator*(const double, const mat&);
	friend mat operator*(const mat&, const double);
	mat operator^(const int) const;
	void operator=(const mat&);

	// friends
	friend vec asVec(const mat&);
	friend unsigned factorial(const unsigned);
	friend mat steppedView(mat);
	friend mat betterSteppedView(const mat&);
	friend mat transpose(const mat&);
	friend int rk(mat);
	friend double det(mat);

	// container
	double** el;

private:	

	unsigned** transpositions();
	int countNonZero(mat);
	void swap(unsigned&, unsigned&);
};

class vec
{
public:

	vec(unsigned);
	vec(const vec&);
	~vec();

	// operators overloading
	double& operator()(const unsigned);
	vec operator+(const vec&) const;
	vec operator-(const vec&) const;
	friend vec operator*(const double, const vec&);
	friend vec operator*(const vec&, const double);
	void operator=(const vec&);
	friend std::ostream& operator<<(std::ostream&, const vec&);

	// turning vec to mat
	friend mat asRow(const vec&);
	friend mat asCol(const vec&);

	// friends
	friend bool areCollinear(const vec&, const vec&);
	friend bool areCoplanar(const vec&, const vec&, const vec&);
	friend double length(const vec&);
	friend double cosBetween(const vec&, const vec&);
	friend double dotProd(const vec&, const vec&);
	friend vec crossProd(const vec&, const vec&);

private:

	mat* m;
};

class complex
{
public:

	complex();
	complex(const double, const double);
	complex(const complex&);
	~complex();

	void operator()(const double, const double);
	complex operator+(const complex&) const;
	complex operator-(const complex&) const;
	complex operator*(const complex&) const;
	complex operator^(const int) const;
	void operator=(const complex&);
	void operator=(const double);
	friend std::ostream& operator<<(std::ostream&, const complex&);

	// friends
	friend double& Re(const complex&);
	friend double& Im(const complex&);
	friend double mod(const complex&);
	friend double arg(const complex&);
	friend complex conj(const complex&);
	friend complex* complRoots(const complex&, const int);

private:

	mat* m;
};

double getInterpolatedValue(double*, double*, unsigned, double);
