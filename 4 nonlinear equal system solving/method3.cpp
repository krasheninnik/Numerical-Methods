#include "stdafx.h"
#include "matrix.h"

double Matrix::derivate(const func &f, realArr x, const double step, const int arg) {
	x[arg] += step;  	
 	double deriv = f(x);
	x[arg] -= 2 * step;
	deriv -= f(x);
	deriv /= 2 * step;
	return deriv;
}


void Matrix::calcSqJacobiMatrix3(const std::vector<bool> &mask, const realArr &x) {
	const double step = 1e-1
		;
	for (int equal = 0, i = 0; equal < numEquals; equal++) {
		if (mask[equal]) {
			for (int j = 0; j < jacobi.size(); j++) {
				jacobi[i][j] = derivate(funcs[equal], x, step, j);
			}
			i++;
		}
	}
}

void Matrix::method3(realArr &x, realArr &currF, realArr &temp,
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
			calcSqJacobiMatrix3(mask, x);

			for (int i = 0, j = 0; i < numEquals; i++) // forming right part temp F:
			{
				if (mask[i]) tempF[j++] = -currF[i];
			}

			methodGaussa(dx, tempF); 					// Solve [Jacobi]*[dx] = [tempF]

			findX2(temp1, x, dx, currF, B, normCurrF);	// add dx to x with B param:


			iter++;
			logInfo(iter, B, normCurrF, x);				// log info
			xPath.push_back(x);
		}
	}
	catch (std::string e) {
		xPath.push_back(x);
		logInfo(++iter, B, normCurrF, x);				// log info
		fout << std::endl << e << std::endl;
		fout.close();
	}
	fout.close();
}
