#pragma once
#include <vector>

class Matrix {
public:
	Matrix();
	~Matrix() {};

	int get_dim();
	void load();

	std::vector<double> method_Jacobi(std::vector<double>&, std::vector<double>&,
		std::vector<double>&, std::vector<double>&, double);
	std::vector<double> method_Gauss_Seidel(std::vector<double>&,
		std::vector<double>&, std::vector<double>&, double);

	std::vector<double> method_block_relaxation(std::vector<double>& x0, 
		std::vector<double>& f, std::vector<double>& temp, double w,
		std::vector<double>&, std::vector<double>&);
	void block_factorization(int, std::vector<double>&, std::vector<double>&);
	std::vector<double> solve_of_block_system(std::vector<double>, std::vector<double>, int,
		std::vector<double>&, std::vector<double>&);


	double calc_sum(int, std::vector<double>&);// sum of multiplicate elems of row matrix with corresponding vector's elems
	double calc_relative_discrepancy(std::vector<double> &, std::vector<double>&, std::vector<double>&);
	std::vector<double> multiplicate_with_vector(std::vector<double>&, std::vector<double>&);



	double estimate_condition_number(std::vector<double> exact_x, std::vector<double> x,
			std::vector<double> &f);
	/*
	// for exploration:1
	
	void explore_matrix(std::vector<double>& x1, std::vector<double> &x0, std::vector<double>& exact_x,
		std::vector<double> &f, std::vector<double> &temp, std::vector<double>& mid_factor_diag1, std::vector<double>& mid_factor_diag2);

	std::vector<double> method_Jacobi_e(std::vector<double>&, std::vector<double>&,
		std::vector<double>&, std::vector<double>&, double);
	std::vector<double> method_Gauss_Seidel_e(std::vector<double>&,
		std::vector<double>&, std::vector<double>&, double);
	std::vector<double> Matrix::method_block_relaxation_e(int, std::vector<double>&,
		std::vector<double>&, std::vector<double>&, double, std::vector<double>&, std::vector<double>&);
		*/

private:
	static const int main = 4;
	static const int diags = 9;

	std::vector<double> low_diag1;
	std::vector<double> low_diag2;
	std::vector<double> low_diag3;

	std::vector<double> mid_diag1;
	std::vector<double> mid_diag2;
	std::vector<double> mid_diag3;

	std::vector<double> high_diag1;
	std::vector<double> high_diag2;
	std::vector<double> high_diag3;
	
	int n, m, block_size, max_iter;
	double accuracy;
};

