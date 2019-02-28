#pragma once
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include <iostream>
#include <assert.h>
#include <set>
#include <string>

#define TWO

void vector_add(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result);
double calc_norm(std::vector<double>& vector);

typedef std::vector<double> realArr;
typedef std::function<double(const realArr&)> func;

class Matrix {
public:
	void load();
	void init();

	int getNumEquals();
	int getNumArgs();
	int getMethod();

	void solve();
	double normInPoint(const realArr&, realArr&);
	std::vector<realArr> xPath;

private:
	void method1(realArr&, realArr&, realArr&, realArr&, realArr&, std::vector<bool>&);
	void method2(realArr&, realArr&, realArr&, realArr&, realArr&, realArr&, std::vector<bool>&);
	void method3(realArr&, realArr&, realArr&, realArr&, realArr&, realArr&, std::vector<bool>&);

	void logInfo(const int, const double, const double, const realArr&);
	void calcCurrF(realArr&, const realArr &);
	void methodGaussa(realArr&, realArr&);
	double derivate(const func &, realArr, const double, const int);

	void calcJacobiMatrix(const realArr &);
	void calcSqJacobiMatrix1(const std::vector<bool> &);
	void calcSqJacobiMatrix2(const std::vector<bool>&, const realArr&);
	void calcSqJacobiMatrix3(const std::vector<bool>&, const realArr&);

	void calcMask1(std::vector<bool>&, realArr&, realArr&, const realArr&);
	void calcMask2(std::vector<bool>&, realArr&, const realArr&);

	void findX1(const std::vector<bool>&,realArr&, realArr&, const realArr&, realArr&, double &, double &);
	void findX2(realArr&, realArr&, const realArr&, realArr&, double&, double&);

private:
	std::ofstream fout;
	std::vector<realArr> elems;
	std::vector<realArr> jacobi;
	
	std::vector<func> funcs;	// functions
	std::vector<std::vector<func>> dfuncs;	// derivative of function

	double eps_b;
	double eps_norm;
	int	max_iter;

	int numEquals;						// amount params
	int numArgs;						// amount equals

	int method;						// solve method 
};