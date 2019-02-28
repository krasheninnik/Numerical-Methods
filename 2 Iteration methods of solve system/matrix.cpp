#include "pch.h"
#include "matrix.h"
#include "stdio.h"
#include <string>
#include <algorithm>
#include <math.h>
#include "conio.h"

double calc_norm(std::vector<double> &x) {
	double norm = 0;
	for (int i = 0; i < x.size(); i++) {
		norm += x[i] * x[i];
	}
	norm = sqrt(norm);
	return norm;
}

std::vector<double> load_vector(int size, std::string filename, std::vector<double> &vect) {
	FILE *in;
	fopen_s(&in, filename.c_str(), "r");
	for (int i = 0; i < size; i++) fscanf_s(in, "%lf ", &vect[i]);
	fclose(in);

	return vect;
}
void out_vector(std::vector<double> &x) {
	FILE *in;
	fopen_s(&in, "x.txt", "w");
	for (int i = 0; i < x.size(); i++) {
		fprintf_s(in, "lf", x[i]);
	}
	fclose(in);
}
Matrix::Matrix()									// initialization
{
	low_diag1 = std::vector<double>(0);
	low_diag2 = std::vector<double>(0);
	low_diag3 = std::vector<double>(0);

	mid_diag1 = std::vector<double>(0);
	mid_diag2 = std::vector<double>(0);
	mid_diag3 = std::vector<double>(0);

	high_diag1 = std::vector<double>(0);
	high_diag2 = std::vector<double>(0);
	high_diag3 = std::vector<double>(0);
}

int Matrix::get_dim() {
	return n; 
}
void Matrix::load() {
	FILE *in;
	fopen_s(&in, "matrix.txt", "r");
	fscanf_s(in, "%d %d %d", &n, &m, &block_size);

	low_diag1.resize(n - m - 4);
	low_diag2.resize(n - m - 3);
	low_diag3.resize(n - m - 2);
	mid_diag1.resize(n - 1);
	mid_diag2.resize(n);
	mid_diag3.resize(n - 1);
	high_diag1.resize(n - m - 2);
	high_diag2.resize(n - m - 3);
	high_diag3.resize(n - m - 4);

	int i;
	for (i = 0; i < low_diag1.size(); i++) fscanf_s(in, "%lf", &low_diag1[i]);
	for (i = 0; i < low_diag2.size(); i++) fscanf_s(in, "%lf", &low_diag2[i]);
	for (i = 0; i < low_diag3.size(); i++) fscanf_s(in, "%lf", &low_diag3[i]);

	for (i = 0; i < mid_diag1.size(); i++) fscanf_s(in, "%lf", &mid_diag1[i]);
	for (i = 0; i < mid_diag2.size(); i++) fscanf_s(in, "%lf", &mid_diag2[i]);
	for (i = 0; i < mid_diag3.size(); i++) fscanf_s(in, "%lf", &mid_diag3[i]);

	for (i = 0; i < high_diag1.size(); i++) fscanf_s(in, "%lf", &high_diag1[i]);
	for (i = 0; i < high_diag2.size(); i++) fscanf_s(in, "%lf", &high_diag2[i]);
	for (i = 0; i < high_diag3.size(); i++) fscanf_s(in, "%lf", &high_diag3[i]);

	fclose(in);

	fopen_s(&in, "accuracy.txt", "r");				// info about accuracy
	fscanf_s(in, "%lf %d", &accuracy, &max_iter);
	fclose(in);
}

double Matrix::calc_sum(int row, std::vector<double> &x) {
	int ldi = row - m - 2;  // first low_diags_index
	int hdi = row + m + 2;	// first high diags_index
	int i;
	double sum = 0;

	// Throw low diags
	if (ldi > 1) {  // use all 3 low diags.
		sum += low_diag1[ldi - 2] * x[ldi - 2];	
		sum += low_diag2[ldi - 1] * x[ldi - 1];
		sum += low_diag3[ldi] * x[ldi];
	} 
	else {
		if (ldi > 0) { // use 2 and 3 low diags.
				sum += low_diag2[ldi - 1] * x[ldi - 1];
				sum += low_diag3[ldi] * x[ldi];
		}
	
		else {
			if (ldi == 0) { // use only 3-th low diags.
				sum += low_diag3[ldi] * x[ldi];
			}
		}
	}

	// Throw main diags
	if (row > 0) {
		if (row < n - 1) { // use all main diags
			sum += mid_diag1[row - 1] * x[row - 1];		
			sum += mid_diag2[row] * x[row];				
			sum += mid_diag3[row] * x[row + 1];			
		}

		else {		// use 1,2 main diags
			sum += mid_diag1[row - 1] * x[row - 1];		
			sum += mid_diag2[row] * x[row];		
		}
	}
	else {			// use 2,3 main diags
		sum += mid_diag2[row] * x[row];	
		sum += mid_diag3[row] * x[row + 1];	
	}

	// Throw high diags 
	if (n - hdi > 2) {	// use all high diags
		sum += high_diag1[row] * x[hdi++];
		sum += high_diag2[row] * x[hdi++];
		sum += high_diag3[row] * x[hdi];
	}
	else {
		if (n - hdi > 1) {	// use 1,2-th high diags
			sum += high_diag1[row] * x[hdi++];
			sum += high_diag2[row] * x[hdi];
		}
		else {			// use only 1-th high diags1
			if (n - hdi > 0) {
				sum += high_diag1[row] * x[hdi];
			}
		}
	}

	return sum;
}

std::vector<double> Matrix::method_Jacobi(std::vector<double> &x1, std::vector<double> &x0,
								std::vector<double> &f, std::vector<double> &temp, double w) {
	printf_s("Method Jacobi: ");
	double relative_discrepancy;

	for (int i = 1; i <= max_iter; i++) {

		for(int j = 0; j < n; j++) {		// Calculate Xi of next iteration
			x1[j] = x0[j] + (w/mid_diag2[j]) * (f[j] - calc_sum(j, x0));
		}

		// check discrepancy
		relative_discrepancy = calc_relative_discrepancy(x1, f, temp);
		printf_s("iteration: %d \t relative discrepancy: %e \n", i, relative_discrepancy);

		if (relative_discrepancy < accuracy) {
			return x1;
		}
		std::swap(x1, x0);		// preparate next iteration			
	}

	printf("i > max_iters!");
	
	return x0;		// return x? or error , amount iter grether then max_iter

}

std::vector<double> Matrix::method_Gauss_Seidel(std::vector<double> &x0,
	std::vector<double> &f, std::vector<double> &temp, double w) {
	printf_s("Method Gauss-Seidel:");
	double relative_discrepancy;
	double temp_sum;

	for (int i = 1; i <= max_iter; i++) {

		for (int j = 0; j < n; j++) {		// Calculate Xi of next iteration
			x0[j] = x0[j] + (w / mid_diag2[j]) * (f[j] - calc_sum(j, x0));
		}

		// check discrepancy
		relative_discrepancy = calc_relative_discrepancy(x0, f, temp);
		printf_s("iteration: %d \t relative discrepancy: %e \n", i, relative_discrepancy);

		if (relative_discrepancy < accuracy) {
		
			return x0;
		}
	}

	printf("i > max_iters!");

	return x0;		// return x? or error , amount iter grether then max_iter

}

double Matrix::calc_relative_discrepancy(std::vector<double> &x, std::vector<double> &f, std::vector<double> &ax) {
	ax = multiplicate_with_vector(x, ax);	    	      // Calculate Ax.
	for (int i = 0; i < ax.size(); i++) ax[i] -= f[i];    // Calculate ( Ax - f )
	double relative_discr = calc_norm(ax) / calc_norm(f); // Calculate discrepancy

	return relative_discr;
}

std::vector<double> Matrix::multiplicate_with_vector(std::vector<double> &x, std::vector<double> &f) {
	for (int i = 0; i < n; i++) f[i] = calc_sum(i, x);
	return f;
}

double  Matrix::estimate_condition_number(std::vector<double> exact_x, std::vector<double> x,
	std::vector<double> &f) {

	double condition_number;
	double norm_exact_x = calc_norm(exact_x);
	double norm_f = calc_norm(f);

	std::vector<double> *discrepancy = &exact_x;

	(*discrepancy) = multiplicate_with_vector(x, (*discrepancy));

	for (int i = 0; i < x.size(); i++) {
		(*discrepancy)[i] -= f[i];
		x[i] -= i + 1;
	}

	condition_number = (calc_norm(x) / norm_exact_x) / (calc_norm((*discrepancy)) / norm_f);
	return condition_number;
}
