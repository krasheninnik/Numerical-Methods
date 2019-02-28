#include "stdafx.h"
#include "matrix.h"

void Matrix::calcSqJacobiMatrix2(const std::vector<bool> &mask, const realArr &x) {
	for (int equal = 0, i = 0; equal < numEquals; equal++) {
		if (mask[equal]) {
			for (int j = 0; j < jacobi.size(); j++) {
				jacobi[i][j] = dfuncs[equal][j](x);
			}
			i++;
		}
	}
}

void Matrix::calcMask2(std::vector<bool> &mask, realArr &temp, const realArr &currF) {

	for (int j = 0; j < mask.size(); j++) mask[j] = 1;	// init mask

	temp = currF;
	std::for_each(temp.begin(), temp.end(), [](double &el) { el = abs(el); });
	std::sort(temp.begin(), temp.end());

	int size_diff = numEquals - numArgs;

	// mask[i] = 0, if equals[i] willn't use in Jacobi matrix 
	for (int i = 0; i < size_diff; i++) {
		for (int j = 0; j < currF.size(); j++) {
			if (temp[i] == abs(currF[j]) && mask[j]) {
				mask[j] = 0;
				break;
			}
		}
	}
}

void Matrix::findX2(realArr &temp, realArr &x,
	const realArr &dx, realArr &currF, double &B, double &normCurrF) {

	double normTempF;
	temp = x;

	B = 2;
	do {
		B /= 2;
		for (int i = 0; i < x.size(); i++) {
			x[i] = temp[i] + B * dx[i];
		}

		calcCurrF(currF, x);
		normTempF = calc_norm(currF);
	} while (normTempF > normCurrF && B > eps_b);

	if (B < eps_b) {
		throw std::string("B < eps_b\n");
	}

	normCurrF = normTempF;
}

void Matrix::method2(realArr &x, realArr &currF, realArr &temp,
	realArr &tempF, realArr &dx, realArr &temp1, std::vector<bool> &mask) {

	fout.open("result.txt");

	assert(mask.size() == numEquals);
	assert(temp.size() == numEquals);
	assert(currF.size() == numEquals);
	assert(temp1.size() == numArgs);
	assert(x.size() == numArgs);
	assert(dx.size() == numArgs);
	assert(tempF.size() == numArgs);

	double B = 0;

	calcCurrF(currF, x);		// calculate current F
	double normF = calc_norm(currF);
	double normCurrF = normF;

	xPath.push_back(x);
	int iter = 0;
	try {
		while (normCurrF / normF > eps_norm && iter < max_iter) {
			calcMask2(mask, temp, currF);
			calcSqJacobiMatrix2(mask, x);

			for (int i = 0, j = 0; i < numEquals; i++) // forming right part temp F:
				{ if (mask[i]) tempF[j++] = -currF[i]; }

			methodGaussa(dx, tempF); 					// Solve [Jacobi]*[dx] = [tempF]

			findX2(temp1, x, dx, currF, B, normCurrF);	// add dx to x with B param:


			iter++;
			logInfo(iter, B, normCurrF, x);				// log info
			xPath.push_back(x);				
		}
	}
	catch(std::string e) {
		xPath.push_back(x);
		logInfo(++iter, B, normCurrF, x);				// log info
		fout << std::endl << e << std::endl;
		fout.close();
	}
	fout.close();
}
