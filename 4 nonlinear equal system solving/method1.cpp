#include "stdafx.h"
#include "matrix.h"


void Matrix::calcJacobiMatrix(const realArr &x) {
	// storage full Jacobi matrix in ELEMS matirx.
	for (int i = 0; i < numEquals; i++) {
		for (int j = 0; j < numArgs; j++) {
			elems[i][j] = dfuncs[i][j](x);
		}
	}
}

void Matrix::calcSqJacobiMatrix1(const std::vector<bool> &mask) {
	for (int arg = 0, j = 0; arg < numArgs; arg++) {
		if (mask[arg]) {
			for (int i = 0; i < jacobi.size(); i++) {
				jacobi[i][j] = elems[i][arg];
			}
			j++;
		}
	}
}

void Matrix::calcMask1(std::vector<bool> &mask, realArr &temp, realArr &maxElems, const realArr &x) {

	for (int j = 0; j < mask.size(); j++) mask[j] = 1;	// init mask
	calcJacobiMatrix(x);

	// find max in columns
	for (int j = 0; j < numArgs; j++) {
		maxElems[j] = abs(elems[0][j]);
		for (int i = 0; i < numEquals; i++) {
			if (abs(elems[i][j]) > maxElems[j]) maxElems[j] = abs(elems[i][j]);
		}
	}

	temp = maxElems;
	std::sort(temp.begin(), temp.end());

	int size_diff = numArgs - numEquals;

	// mask[i] = 0, if args[i] willn't use in squad Jacobi matrix 
	for (int i = 0; i < size_diff; i++) {
		for (int j = 0; j < maxElems.size(); j++) {
			if (temp[i] == maxElems[j] && mask[j]) {
				mask[j] = 0;
				break;
			}
		}
	}
}

void Matrix::findX1(const std::vector<bool> &mask, realArr &temp, realArr &x, 
	const realArr &dx, realArr &currF, double &B, double &normCurrF) {

	double normTempF;
	temp = x;

	B = 2;
	do {
		B /= 2;
		for (int i = 0, j = 0; i < x.size(); i++) {
			if(mask[i]) x[i] = temp[i] + B * dx[j++];
		}

		calcCurrF(currF, x);
		normTempF = calc_norm(currF);

	} while (normTempF > normCurrF && B > eps_b);

	if (B < eps_b) {
		throw std::string("B < eps_b\n");
	}

	normCurrF = normTempF;
}

void Matrix::method1(realArr &x, realArr &currF, realArr &temp,
	realArr &maxElems, realArr &dx, std::vector<bool> &mask) {

	fout.open("result.txt");

	assert(x.size() == numArgs);
	assert(mask.size() == numArgs);
	assert(temp.size() == numArgs);
	assert(maxElems.size() == numArgs);

	assert(currF.size() == numEquals);
	assert(dx.size() == numEquals);
	
	double B = 0;

	calcCurrF(currF, x);		// calculate current F
	double normF = calc_norm(currF);
	double normCurrF = normF;

	xPath.push_back(x);
	int iter = 0;
	try {
		while (normCurrF / normF > eps_norm && iter < max_iter) {
			calcMask1(mask, temp, maxElems, x);		// calc mask
			calcSqJacobiMatrix1(mask);				// form Jacobi matrix:

			for (int i = 0; i < currF.size(); i++) currF[i] *= -1;	// form Right part
			methodGaussa(dx, currF);				// Solve [Jacobi]*[dx] = [currF]

			findX1(mask, temp, x, dx, currF, B, normCurrF);	// add dx to x with B param:

			iter++;
			logInfo(iter, B, normCurrF, x);				// log info
			xPath.push_back(x);
		}
	}
	catch (std::string e) {
		xPath.push_back(x);
		logInfo(iter, B, normCurrF, x);				// log info
		std:fout << std::endl << e << std::endl;
		fout.close();
		////for (auto &x : xPath) {
		//	for (auto &coord : x) fout << coord;
		//	fout << '\n';
		//}
	}
	fout.close();
}