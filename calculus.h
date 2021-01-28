#pragma once
#include <iostream>
#include <cmath>

const double PI = 3.141592653589793;
bool comp_d(const double, const double);

class mat
{
public:

	mat(const unsigned, const unsigned);
	mat(const mat&);
	~mat();

	unsigned rows, cols;

	double& operator()(const unsigned, const unsigned);
	mat operator+(const mat&) const;
	mat operator-(const mat&) const;
	mat operator*(const mat&) const;
	friend mat operator*(const double, const mat&);
	friend mat operator*(const mat&, const double);
	mat operator^(const int) const;
	void operator=(const mat&);

	// friends
	friend unsigned factorial(const unsigned);
	friend mat steppedView(mat);
	friend mat betterSteppedView(const mat&);
	friend mat transpose(const mat&);
	friend int rk(mat);
	friend double det(mat);

private:

	double** el;

	int countNonZero(mat);
	void swap(unsigned&, unsigned&);
	unsigned** transpositions();
	friend mat roundByEps(const mat&);
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

	// friends
	friend std::ostream& operator<<(std::ostream&, const complex&);
	friend double& Re(const complex&);
	friend double& Im(const complex&);
	friend double mod(const complex&);
	friend double arg(const complex&);
	friend complex conj(const complex&);
	friend complex* complRoots(const complex&, const int);

private:

	mat* m;
};
