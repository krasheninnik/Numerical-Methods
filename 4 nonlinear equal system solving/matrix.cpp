#include "stdafx.h"
#include "matrix.h"

void save_vector(std::vector<double> &vect, std::string filename, std::string frm);
void load_vector(int size, std::string filename, std::vector<double> &vect);
double calc_norm(std::vector<double>& vector);

int Matrix::getNumEquals() { return numEquals; }
int Matrix::getNumArgs() { return numArgs; }
int Matrix::getMethod() { return method; }

void Matrix::logInfo(const int iter, const double B, const double normF, const realArr &x) {
	fout.precision(11);
	fout <<	"iter: " << iter << " B: " << B << " norm: " << normF << "\t x: ";

	for (auto coord : x) fout << coord << "  ";
	fout << std::endl;
}

void Matrix::load() {
	std::ifstream fin("accuracy.txt");
	fin >> eps_norm >> eps_b >> max_iter;
	fin.close();
		
	
	fin.open("matrix.txt");
	fin >> numEquals >> numArgs;
	fin.close();

	fin.open("method.txt");
	fin >> method;
	fin.close();

	if (method < 1 || method > 3)
		std::cout << "err: indefinite solve method" << std::endl;
}

void Matrix::calcCurrF(realArr &currF, const realArr &x) {
	for (int i = 0; i < numEquals; i++) {
		currF[i] = funcs[i](x);
	}
}

void Matrix::methodGaussa(realArr& x, realArr& f) {
	// row reducing
	int maxElInd;
	double maxEl;
	double coef;

	for (int i = 0; i < jacobi.size() - 1; i++) {
		// select main elems from column:
		maxEl = abs(jacobi[i][i]);
		maxElInd = i;

		for (int j = i + 1; j < jacobi.size(); j++) {
			if (abs(jacobi[j][i]) > maxEl) {
				maxEl = abs(jacobi[j][i]);
				maxElInd = j;
			}
		}

		if (maxEl == 0) {
			throw std::string("Can't solve. Zero on diag");
		}

		// row swap
		if (maxElInd != i) {
			std::swap(jacobi[i], jacobi[maxElInd]);
			std::swap(f[i], f[maxElInd]);
		}

		// zeroing column
		for (int j = i + 1; j < jacobi.size(); j++) {
			if(jacobi[j][i] != 0) {
				coef = jacobi[j][i] / jacobi[i][i];

				for (int k = i; k < jacobi.size(); k++) {
					jacobi[j][k] -= coef * jacobi[i][k];
				}

				f[j] -= coef * f[i];
			}
		}
	}

	if(jacobi[jacobi.size() - 1][jacobi.size() - 1] == 0) 
		throw std::string("Can't solve. Zero on diag");

	// backward Gauss:

	double sum = 0;
	for (int i = jacobi.size() - 1; i >= 0; i--) {
		sum = 0;
		for (int j = jacobi.size() - 1; j > i; j--) {
			sum += jacobi[i][j] * x[j];
		}

		x[i] = (f[i] - sum) / jacobi[i][i];
	}
}
void Matrix::solve() {
	realArr x0(getNumArgs());
	load_vector(getNumArgs(), "x0.txt", x0);

	int sq_size = min(getNumEquals(), getNumArgs());

	// memory allocate for methods
	realArr resourceEq1(getNumEquals()), resourceEq2(getNumEquals()),
		resourceArg1(getNumArgs()), resourceArg2(getNumArgs()), resourceArg3(getNumArgs());

	std::vector<bool> mask(getMethod() == 1 ? getNumArgs() : getNumEquals());

	switch (method) {
	case 1:		// equals < args 
		method1(x0, resourceEq1, resourceArg1, resourceArg2, resourceEq2, mask);
		break;

	case 2:		// equals >= args 
		method2(x0, resourceEq1, resourceEq2, resourceArg1, resourceArg2, resourceArg3, mask);
		break;

	case 3:		// equals >= args; numeric calc derivative.
		method3(x0, resourceEq1, resourceEq2, resourceArg1, resourceArg2, resourceArg3, mask);
		break;
	}

	save_vector(x0, "1x.txt", "%e");
}

double Matrix::normInPoint(const realArr& point, realArr& f) {
	calcCurrF(f, point);
	return calc_norm(f);
}