#pragma once
#include "pch.h"
#include <vector>

typedef std::vector<double> realArr;

//#define ONE
class Spline {			// Cubic Hermite spline
public:
	void load();
	void calcDerivates();
	void calcSpline();
	void output();

private:
	double n;
	realArr x;			// point x
	realArr f;			// value in point x
	realArr df;			// derivate in point 

	realArr spline;		// calculated Spline values

	int partitCoef;		// amount calculated values in each [xi, xi+1]
};
