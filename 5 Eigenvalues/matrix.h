#pragma once
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include <iostream>
#include <assert.h>
#include <set>
#include <string>
#include <stdint.h>

double calc_norm(std::vector<double>& vector);

typedef std::vector<double> realArr;

class Matrix {
public:
	void load();
	void findEigens();

private:
	void LUdecompose();
	void solve(realArr&, realArr&, realArr&);
	void multiplicate(realArr&, realArr&);

	void  vectorNormalization(realArr&);
	bool  checkFinding(realArr&, realArr&, double&);

	double findMaxEigenValue(realArr&, realArr&, uint64_t &);
	double findMinEigenValue(realArr&, realArr&, realArr&, uint64_t &);
	
private:
	unsigned n;		 // size of matrix
	std::vector<realArr> elems;

	double eps;		// error of eighValue
	double beta;	// error of close to zero numbers
};