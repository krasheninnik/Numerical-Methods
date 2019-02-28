#include "pch.h"
#include "spline.h"

void Spline::load() {
	/* input data for spline and memory allocation*/
	std::ifstream fin("data.txt");
	double point, value;
	while (fin >> point >> value) {
		x.push_back(point);
		f.push_back(value);		
	}
	fin.close();

	fin.open("accuracy.txt");
	if (fin.is_open()) fin >> partitCoef;
	else partitCoef = 100;
	
	/* memory allocation for spline */
	df.resize(f.size());
	spline.resize((x.size() - 1) * partitCoef + 1);
}

void Spline::output() {
	std::ofstream fout("spline.txt");
	//std::ofstream fout(R"(C:\Users\krasheninnik.2016\PycharmProjects\plotSpline\spline.txt)");
	fout.precision(16);

	double hk;		// length of [xi, xi+1];
	double step;	// step to calculate value of spline in [xi, xi+1]
	double _x;		// point to calculate value of spline

	int j;			// index of spline calculated elem

	for (int i = 0; i < x.size() - 1; i++) {
		hk = x[i + 1] - x[i];
		step = hk / partitCoef;

		j = i * partitCoef;	// index of next spline calculated elem	
		for (int k = 0; k < partitCoef; k++) {
			_x = x[i] + step * k;
			fout << std::setw (20) << _x << "  " << std::setw(20) << spline[j++] << std::endl;
			//fout << _x << ","  << spline[j++] << std::endl;
		}
	}

	// last slpline elem
	_x = x[x.size() - 1];
	fout << std::setw(20) << _x << "  " << std::setw(20) << spline[j++] << std::endl;
	//fout << _x << "  " << spline[j++] << std::endl;

}

inline double derivate(const double& xLeft, const double& xRight, const double& fLeft, const double& fRight) {
	return (fRight - fLeft) / (xRight - xLeft);
}


void Spline::calcDerivates() {
	df[0] = derivate(x[0], x[1], f[0], f[1]);

	for (unsigned i = 1; i < df.size() - 1; i++) {
		df[i] = derivate(x[i - 1], x[i + 1], f[i - 1], f[i + 1]);
	}

	auto lastPos = df.size() - 1;
	df[lastPos] = derivate(x[lastPos - 1], x[lastPos], f[lastPos - 1], f[lastPos]);
}

inline double toLocalCoord(const double&x, const double&xk, const double hk) {
	return (x - xk)/hk;
}

inline double basic1(const double& ksi) {
	return 1 - 3 * pow(ksi, 2) + 2 * pow(ksi, 3);
}

inline double basic2(const double& ksi) {
	return (ksi - 2 * pow(ksi, 2) + pow(ksi, 3));
}

inline double basic3(const double& ksi) {
	return 3 * pow(ksi, 2) - 2 * pow(ksi, 3);
}

inline double basic4(const double& ksi) {
	return -pow(ksi, 2) + pow(ksi, 3);
}


void Spline::calcSpline() {
	calcDerivates(); // prepare calculation of Spline.

	double hk;		// length of [xi, xi+1];
	double step;	// step to calculate value of spline in [xi, xi+1]
	double _x;		// point to calculate value of spline
	double ksi;		// local coordinate

	int j = 0;			// index of spline calculated elem

	for (int i = 0; i < x.size() - 1; i++) {
		hk = x[i + 1] - x[i];
		step = hk / partitCoef;

		for (int k = 0; k < partitCoef; k++) {
			_x = x[i] + step * k;

			ksi = toLocalCoord(_x, x[i], hk);
				
			spline[j++] = f[i]			 * basic1(ksi) +
						  df[i]		* hk * basic2(ksi) +
						  f[i + 1]		 * basic3(ksi) +
						  df[i + 1] * hk * basic4(ksi);
		}
	}

	// last spline value:
	auto lastPos = x.size() - 1;
	ksi = 1;
	_x = x[lastPos];
	spline[j] =  f[lastPos - 1]		 * basic1(ksi) +
				df[lastPos - 1] * hk * basic2(ksi) +
				 f[lastPos]			 * basic3(ksi) +
				df[lastPos]		* hk * basic4(ksi);

}