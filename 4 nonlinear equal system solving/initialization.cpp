#include "stdafx.h"
#include "matrix.h"


double  q = 0.000001;	// weight for some equals

void Matrix::init() {

#ifdef ONE
	// две непересекающиеся окружности
	method = 2;
	numEquals = 2;
	numArgs = 2;

	//set functions //
	funcs = std::vector<func>(numEquals);

	funcs[0] = [](const realArr &x) {return pow(x[0], 2) + pow(x[1], 2) - 1; };
	funcs[1] = [](const realArr &x) {return pow(x[0] - 3, 2) + pow(x[1] - 2, 2) - 4; };

	// set partial derivatieves of functions // 
	dfuncs = std::vector<std::vector<func>>(numEquals);
	for (auto &df : dfuncs) df = std::vector<func>(numArgs);

	dfuncs[0][0] = [](const realArr &x) {return 2 * x[0]; };
	dfuncs[0][1] = [](const realArr &x) {return 2 * x[1]; };

	dfuncs[1][0] = [](const realArr &x) {return 2 * (x[0] - 3); };
	dfuncs[1][1] = [](const realArr &x) {return 2 * (x[1] - 2); };
	
#endif

#ifdef TWO
	// две непересекающиеся окружности
	method = 3;
	numEquals = 2;
	numArgs = 2;

	//set functions //
	funcs = std::vector<func>(numEquals);

	funcs[0] = [](const realArr &x) {return pow(x[0] + 2, 2) + pow(x[1] + 2, 2) - 4; };
	funcs[1] = [](const realArr &x) {return pow(x[0] - 1, 2) + pow(x[1] - 2, 2) - 9; };

	// set partial derivatieves of functions // 
	dfuncs = std::vector<std::vector<func>>(numEquals);
	for (auto &df : dfuncs) df = std::vector<func>(numArgs);

	dfuncs[0][0] = [](const realArr &x) {return 2 * (x[0] + 2); };
	dfuncs[0][1] = [](const realArr &x) {return 2 * (x[1] + 2); };

	dfuncs[1][0] = [](const realArr &x) {return 2 * (x[0] - 1); };
	dfuncs[1][1] = [](const realArr &x) {return 2 * (x[1] - 2); };
#endif

#ifdef THREE
	// две непересекающиеся окружности
	method = 3;
	numEquals = 2;
	numArgs = 2;

	//set functions //
	funcs = std::vector<func>(numEquals);

	funcs[0] = [](const realArr &x) {return pow(x[0] + 2, 2) + pow(x[1] + 2, 2) - 9; };
	funcs[1] = [](const realArr &x) {return pow(x[0] - 1, 2) + pow(x[1] - 2, 2) - 9; };

	// set partial derivatieves of functions // 
	dfuncs = std::vector<std::vector<func>>(numEquals);
	for (auto &df : dfuncs) df = std::vector<func>(numArgs);

	dfuncs[0][0] = [](const realArr &x) {return 2 * (x[0] + 2); };
	dfuncs[0][1] = [](const realArr &x) {return 2 * (x[1] + 2); };

	dfuncs[1][0] = [](const realArr &x) {return 2 * (x[0] - 1); };
	dfuncs[1][1] = [](const realArr &x) {return 2 * (x[1] - 2); };
#endif


#ifdef FOUR
	// три непересекающиеся окружности
	method = 2;
	numEquals = 3;
	numArgs = 2;

	//set functions //
	funcs = std::vector<func>(numEquals);

	funcs[0] = [](const realArr &x) {return pow(x[0], 2) + pow(x[1], 2) - 1; };
	funcs[1] = [](const realArr &x) {return pow(x[0] - 3, 2) + pow(x[1] - 2, 2) - 4; };
	funcs[2] = [](const realArr &x) {return pow(x[0] - 3, 2) + pow(x[1] + 3, 2) - 2.25; };

	// set partial derivatieves of functions // 
	dfuncs = std::vector<std::vector<func>>(numEquals);
	for (auto &df : dfuncs) df = std::vector<func>(numArgs);

	dfuncs[0][0] = [](const realArr &x) {return 2 * x[0]; };
	dfuncs[0][1] = [](const realArr &x) {return 2 * x[1]; };

	dfuncs[1][0] = [](const realArr &x) {return 2 * (x[0] - 3); };
	dfuncs[1][1] = [](const realArr &x) {return 2 * (x[1] - 2); };

	dfuncs[2][0] = [](const realArr &x) {return 2 * (x[0] - 3); };
	dfuncs[2][1] = [](const realArr &x) {return 2 * (x[1] + 3); };

#endif

#ifdef FIVE
	// три непересекающиеся окружности
	// method = 2;
	// numEquals = 3;
	// numArgs = 2;

	//set functions //
	funcs = std::vector<func>(numEquals);

	//double  q = 15;

	funcs[0] = [](const realArr &x) {return (pow(x[0], 2) + pow(x[1], 2) - 1); };
	funcs[1] = [](const realArr &x) {return pow(x[0] - 3, 2) + pow(x[1] - 2, 2) - 4; };
	funcs[2] = [](const realArr &x) {return q*(pow(x[0] - 3, 2) + pow(x[1] + 3, 2) - 2.25); };

	// set partial derivatieves of functions // 
	dfuncs = std::vector<std::vector<func>>(numEquals);
	for (auto &df : dfuncs) df = std::vector<func>(numArgs);

	dfuncs[0][0] = [](const realArr &x) {return (2 * x[0]); };
	dfuncs[0][1] = [](const realArr &x) {return (2 * x[1]); };

	dfuncs[1][0] = [](const realArr &x) {return 2 * (x[0] - 3); };
	dfuncs[1][1] = [](const realArr &x) {return 2 * (x[1] - 2); };

	dfuncs[2][0] = [](const realArr &x) {return q*(2 * (x[0] - 3)); };
	dfuncs[2][1] = [](const realArr &x) {return q*(2 * (x[1] + 3)); };

#endif

#ifdef SIX
	method = 2;
	numEquals = 2;
	numArgs = 2;

	// set functions //
	funcs = std::vector<func>(numEquals);

	//double  q = 15;

	funcs[0] = [](const realArr &x) {return (pow(x[0], 2) + pow(x[1], 2) - 4); };
	funcs[1] = [](const realArr &x) {return (pow(x[0] - 4, 2) + pow(x[1], 2) - 4); };
	

	// set partial derivatieves of functions // 
	dfuncs = std::vector<std::vector<func>>(numEquals);
	for (auto &df : dfuncs) df = std::vector<func>(numArgs);

	dfuncs[0][0] = [](const realArr &x) {return (2 * x[0]); };
	dfuncs[0][1] = [](const realArr &x) {return (2 * x[1]); };

	dfuncs[1][0] = [](const realArr &x) {return (2 * (x[0] - 4)); };
	dfuncs[1][1] = [](const realArr &x) {return (2 * x[1]); };


#endif
	/*
	//TEST METHOD 2 
	//about matrix 
	method = 2;

	numEquals = 3;
	numArgs   = 2;

	//set functions //
	funcs = std::vector<func>(numEquals);
	// x[0] = x
	// x[1] = y
	funcs[0] = [](const realArr &x) {return x[0]*x[0] + x[1]*x[1] - 4; };			// x^2 + y^2 - 2 = 0
	funcs[1] = [](const realArr &x) {return x[1] - 2; };							// y - 4 = 0
	funcs[2] = [](const realArr &x) {return x[1] - 0.5 *(x[0] - 2)*(x[0] - 2); };	// y - 1/2(x-2)^2

	// set partial derivatieves of functions // 
	dfuncs = std::vector<std::vector<func>>(numEquals);
	for (auto &df : dfuncs) df = std::vector<func>(numArgs);

	dfuncs[0][0] = [](const realArr &x) {return 2*x[0]; };
	dfuncs[0][1] = [](const realArr &x) {return 2*x[1]; };

	dfuncs[1][0] = [](const realArr &x) {return 0; };
	dfuncs[1][1] = [](const realArr &x) {return 1; };

	dfuncs[2][0] = [](const realArr &x) {return -x[0] + 2; };
	dfuncs[2][1] = [](const realArr &x) {return 1; };
	*/

#ifdef SEVEN
	// Test method 1;
	// about matrix
	method = 1;

	numEquals = 2;
	numArgs = 3;

	//set functions //
	funcs = std::vector<func>(numEquals);
	// x[0] = x
	// x[1] = y
	funcs[0] = [](const realArr &x) {return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 4; };	// сфера радиуса 2
	funcs[1] = [](const realArr &x) {return      x[2] - x[0]*x[0] - x[1]*x[1] - 2; };		// параболоид.верш(0.0.2)

	// set partial derivatieves of functions // 
	dfuncs = std::vector<std::vector<func>>(numEquals);
	for (auto &df : dfuncs) df = std::vector<func>(numArgs);

	dfuncs[0][0] = [](const realArr &x) {return 2 * x[0]; };
	dfuncs[0][1] = [](const realArr &x) {return 2 * x[1]; };
	dfuncs[0][2] = [](const realArr &x) {return 2 * x[2]; };

	dfuncs[1][0] = [](const realArr &x) {return -2 * x[0]; };
	dfuncs[1][1] = [](const realArr &x) {return -2 * x[1]; };
	dfuncs[1][2] = [](const realArr &x) {return 1; };
#endif

	/* memory allocation for elems */
	elems = std::vector<realArr>(numEquals);
	for (auto &eq : elems)  eq = std::vector<double>(numArgs);

	int sq_size = min(numEquals, numArgs);
	jacobi = std::vector<realArr>(sq_size);
	for (auto &eq : jacobi) eq = std::vector<double>(sq_size);
}
