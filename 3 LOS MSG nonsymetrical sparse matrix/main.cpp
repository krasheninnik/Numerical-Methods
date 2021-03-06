// CM_3.cpp : Defines the entry point for the console application.
//
#include "pch.h"
#include "not_symmetrical_sparse_matrix.h"
#include "windows.h"
#include <time.h>


void load_vector(int size, std::string filename, std::vector<double> &vect);
void save_vector(std::vector<double> &vect, std::string filename, std::string frm);

int main()
{
	not_symmetrical_sparse_matrix M;
	M.load();

	std::vector<double> x0(M.get_dim()), f(M.get_dim());
	load_vector(f.size(), "pr.txt", f);
	load_vector(f.size(), "x0.txt", x0);

	std::vector<double> r(M.get_dim()),
		z(M.get_dim()),
		p(M.get_dim()),
		temp(M.get_dim()),
		temp2(M.get_dim()),
		ud(M.get_dim()),
		l(M.get_amount_elems()),
		u(M.get_amount_elems());

#ifdef CALCULATE
	// in MCG.. 'p' is temp storage for 'x', which have min discrepancy.
	switch (M.get_method()) {
		case 1:	M.LOS(x0, f, r, z, p, temp); break;
		case 2: M.LOS_LU(l, ud, u, x0, f, r, z, p, temp, temp2); break;
		case 3: M.LOS_DIAG(ud, x0, f, r, z, p, temp, temp2); break;
		case 4: M.MCG(x0, p, f, r, z, temp, temp2); break;
		case 5: M.MCG_LU(l, ud, u, x0, p, f, r, z, temp, temp2); break;
		case 6: M.MCG_DIAG(ud, x0, p, f, r, z, temp, temp2); break;
		default: printf_s("method error\n");
	}

	save_vector(x0, "1result.txt", "%.15e\n");
	for (int i = 0; i < z.size(); i++) z[i] = i + 1 - x0[i];
	save_vector(z, "1error.txt", "%.15e\n");
#else
	const int ITERS = 5;
	not_symmetrical_sparse_matrix test;
	std::vector<double> source_x0(x0);
	std::vector<double> source_di(ud);
	std::vector<double> source_factorizate(l);

	double CalcTimeStart = clock() / (double)CLOCKS_PER_SEC;

	for (int i = 0; i < ITERS; i++) {
		// in MCG.. 'p' is temp storage for 'x', which have min discrepancy.
		test.copy(M);
		x0 = source_x0;
		l = source_factorizate;
		u = source_factorizate;
		ud = source_di;

		switch (M.get_method()) {
			case 1:	M.LOS(x0, f, r, z, p, temp); break;
			case 2: M.LOS_LU(l, ud, u, x0, f, r, z, p, temp, temp2); break;
			case 3: M.LOS_DIAG(ud, x0, f, r, z, p, temp, temp2); break;
			case 4: M.MCG(x0, p, f, r, z, temp, temp2); break;
			case 5: M.MCG_LU(l, ud, u, x0, p, f, r, z, temp, temp2); break;
			case 6: M.MCG_DIAG(ud, x0, p, f, r, z, temp, temp2); break;
			default: printf_s("method error\n");
		}
	}

	double CalcTimeStop = clock() / (double)CLOCKS_PER_SEC;

	double CopyTimeStart = clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < ITERS; i++) {
		test.copy(M);
		x0 = source_x0;
		u = source_factorizate;
		l = source_factorizate;
		ud = source_di;
	}
	double CopyTimeStop = clock() / (double)CLOCKS_PER_SEC;

	printf("\nEx Time: %f sec\n", (CalcTimeStop - CalcTimeStart) / ITERS);
	printf("Copy Time: %f sec\n", (CopyTimeStop - CopyTimeStart) / ITERS);
	printf("Time: %f sec\n", ((CalcTimeStop - CalcTimeStart) - (CopyTimeStop - CopyTimeStart))/ITERS);
	test.~not_symmetrical_sparse_matrix();

#endif



	M.~not_symmetrical_sparse_matrix();
    return 0;
}
