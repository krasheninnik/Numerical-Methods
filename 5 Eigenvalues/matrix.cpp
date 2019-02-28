#include "pch.h"
#include "matrix.h"

void Matrix::load() {
	/* input size of matrix*/
	std::ifstream fin("matrix.txt");
	fin >> n;

	/* memory allocation for elems */
	elems = std::vector<realArr>(n);
	for (auto &eq : elems)  eq = std::vector<double>(n);

	/* input elems of matrix */
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = 0; j < n; j++) {
			fin >> elems[i][j];
		}
	}

	fin.close();
	fin.open("accuracy.txt");
	fin >> eps >> beta;
	fin.close();
}

void Matrix::multiplicate(realArr& x, realArr& f) {
	for (unsigned int row = 0; row < n; row++) {
		f[row] = elems[row][0] * x[0];			// initiation f[row] 
		for (unsigned int column = 1; column < n; column++) {
			f[row] += elems[row][column] * x[column];
		}
	}
}

void Matrix::LUdecompose() {
	double sum;

	if (elems[0][0] == 0) throw std::string("Zero on diag. LUdecomposite\n");

	for (unsigned i = 1; i < n; i++) {
		// elems of low trinagle:
		for (unsigned j = 0; j < i; j++) {
			sum = 0;
			for (unsigned k = 0; k < j; k++) sum += elems[i][k] * elems[k][j];

			elems[i][j] -= sum; // elems[j][j];
			elems[i][j] /= elems[j][j];
		}
	
		// elems of high trinagle:
		for (unsigned j = 0; j < i; j++) {
			sum = 0;
			for (unsigned k = 0; k < j; k++) sum += elems[j][k] * elems[k][i];

			elems[j][i] -= sum;
		}

		// diag elems:
		sum = 0;
		for (unsigned k = 0; k < i; k++) sum += elems[i][k] * elems[k][i];
		if (abs(elems[i][i] - sum) < abs(elems[i][i] * 1e-15)) throw std::string("Zero on diag. LUdecomposite\n");
		elems[i][i] -= sum;
	}
}

void Matrix::solve(realArr &x, realArr &f, realArr &b) {
	// LUx = f
	//realArr b; //temp

	// forward steps:	// Lb = f
	double sum;
	for (int i = 0; i < n; i++) {
		sum = 0;
		for (int j = 0; j < i; j++) {
			sum += elems[i][j] * b[j];
		}
		b[i] = f[i] - sum;
	}

	// backward steps:	Ux = b
	for (int i = n - 1; i >= 0; i--) {	
		sum = 0;
		for (int j = i + 1; j < n; j++) {
			sum += elems[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / elems[i][i];
	}
}

inline void Matrix::vectorNormalization(realArr &y) {
	double normY = calc_norm(y);
	for (unsigned i = 0; i < y.size(); i++) y[i] /= normY;
}

void Matrix::findEigens() {
	realArr y(n);	
	realArr temp(n, 1);
	realArr temp1(n, 1);

	for (auto &elem : y) elem = std::rand(); // initializate random y;


	for (auto &elem : y) elem = 1; // initializate random y;


	uint64_t iter = 0;
	temp1 = y;	// save y0;

	try {
		std::cout.precision(16);
		std::cout << "Max Eigen Value: " << findMaxEigenValue(y, temp, iter) << "\titer: "<<  iter <<std::endl;
		std::cout << "Min Eigen Value: " << findMinEigenValue(y = temp1, temp, temp1, iter) 
							<< "\t\titer: " << iter << std::endl;	//recovery y0 (y = temp1)
	}
	catch (std::string e) { std::cout << e;	}
}

bool Matrix::checkFinding(realArr& y0, realArr& y1, double &maxEigh) {
	double tempEigh;
	bool finding = true;
	bool isEigh = false;

	unsigned i = 0;

	for (; i < y0.size() && !isEigh; i++) {  // avoid dividion on zero
		if (abs(y0[i]) < beta) continue;

		isEigh = true;
		maxEigh = y1[i] / y0[i];
	}

	if (isEigh == false) throw std::string("imposible find maxEighvalue with the y0");

	for (; i < y0.size() && finding; i++) {
		if (abs(y0[i]) < beta) continue;

		tempEigh = y1[i] / y0[i];
		if (abs(maxEigh - tempEigh) > abs(maxEigh) * eps) finding = false;
	}

	if (finding) {
		maxEigh = 0;
		unsigned count = 0;

		for (i = 0; i < y0.size(); i++) {
			if (abs(y0[i]) < beta) continue;
			maxEigh += y1[i] / y0[i];
			count++;
		}

		maxEigh /= count;
	}

	return finding;
}

double Matrix::findMaxEigenValue(realArr &y, realArr &temp, uint64_t &iter) {
	double maxEigh, previousEigh = 0;
	double normY, normTemp;

	bool findedMaxEigh = false;
	bool equalRelation = false;
	y[0] = 1;
	y[1] = 1;
	iter = 0;
	while (!findedMaxEigh) {
		vectorNormalization(y);

		multiplicate(y, temp);

		equalRelation = checkFinding(y, temp, maxEigh);


		if (equalRelation) {
			if (abs(maxEigh - previousEigh) < abs(maxEigh) * eps) findedMaxEigh = true;
			previousEigh = maxEigh;
		}

		swap(y, temp);		// prepare next iteration
		iter++;
	}

	return maxEigh;
}

double Matrix::findMinEigenValue(realArr &y, realArr &temp, realArr &temp1, uint64_t &iter) {
	double maxEigh = 0;
	double normY, normTemp;

	bool findedMaxEigh = false;

	LUdecompose();

	bool equalRelation = false;
	double previousEigh = 0;

	iter = 0;
	while (!findedMaxEigh) {

		vectorNormalization(y);

		solve(temp, y, temp1);
		
		equalRelation = checkFinding(y, temp, maxEigh);

		//std::cout << maxEigh << "\r";

		if (equalRelation) {
			if (abs(maxEigh - previousEigh) < abs(maxEigh)*eps) findedMaxEigh = true;
			previousEigh = maxEigh;
		}

		swap(y, temp);		// prepare next iteration
		iter++;
	}

	return 1 / maxEigh;
}