#include "pch.h"
#include "matrix.h"
#include "conio.h"

std::vector<double> load_vector(int size, std::string filename, std::vector<double>& vect);

int main()
{
	Matrix A; 
	A.load();

	int size = A.get_dim();
	std::vector<double> x1(size), x0(size), exact_x(size);
	std::vector<double> f(size), temp(size);
	std::vector<double> mfd1(size - 1), mfd2(size);	// mid_fatorized_diagN

	exact_x = load_vector(size, "x.txt", exact_x);
	f = A.multiplicate_with_vector(exact_x, f);	// calc vecotor F for explore.

	//////////////////////////////////////////////////
	//A.explore_matrix(x1, x0, exact_x, f, temp, mfd1, mfd2);
	/////////////////////////////////////////////////

	//double w = 1.06;
	double w = 1.76;
	
	x1 = A.method_block_relaxation(x0, f, temp, w, mfd1, mfd2);
	//printf("%e", A.estimate_condition_number(exact_x, x1, f));
	
	FILE *out;
	fopen_s(&out, "1result.txt", "w");
	for (int i = 0; i < x1.size(); i++) fprintf(out, "%.15lf\n", x1[i]);
	fclose(out);
	fopen_s(&out, "1error.txt", "w");
	for (int i = 0; i < x1.size(); i++) fprintf(out, "%.8e\n", x1[i] - i - 1);
	fclose(out);

	return 0;
}

