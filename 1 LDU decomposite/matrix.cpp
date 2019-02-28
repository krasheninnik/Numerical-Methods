#pragma once
#include "matrix.h"
#include "stdio.h"
#include "algorithm"
#include "stdio.h"
#include "limits"
#include "conio.h"
#include "cmath"

std::string frmt = "%f ";
std::string frmt_in = "%lf ";

Matrix::Matrix() {};
Matrix::~Matrix() {
	ia.clear();
	au.clear();
	al.clear();
	di.clear();
}

void Matrix::load() {
	// load matrix from files
	int i;
	FILE *in;
	fopen_s(&in, "info.txt", "r");
	fscanf_s(in, "%d %d", &n, &m);
	fclose(in);

	fopen_s(&in, "ia.txt", "r");
	ia.resize(n+1);
	for (i = 0; i < n + 1; i++) fscanf_s(in, "%d", &ia[i]);
	fclose(in);

	fopen_s(&in, "di.txt", "r");
	di.resize(n);
	for (i = 0; i < n; i++) fscanf_s(in, frmt_in.c_str(), &di[i]);
	fclose(in);

	fopen_s(&in, "au.txt", "r");
	au.resize(m);
	for (i = 0; i < m; i++) fscanf_s(in, frmt_in.c_str(), &au[i]);
	fclose(in);

	fopen_s(&in, "al.txt", "r");
	al.resize(m);
	for (i = 0; i < m; i++) fscanf_s(in, frmt_in.c_str(), &al[i]);
	fclose(in);
}

std::vector<mytype> load_vector(int size, std::string filename, std::vector<mytype>& vect) {
	FILE *in;
	fopen_s(&in, filename.c_str(), "r");
	for (int i = 0; i < size; i++) fscanf_s(in, frmt_in.c_str(), &vect[i]);
	fclose(in);

	return vect;
}

int Matrix::get_dim() { return n; }

void Matrix::LDU_decompose() {
	/* 
	 al -> L
	 au -> U
	 di -> di
	*/

	
	double eps;
	eps = typeid(di[0]) == typeid(double) ? 1e-15 : 1e-7;	// calculate eps 

	mytype2 sumL, sumU;
	int i, j, k;
	int row, column, offset, amount_elems;
	int index_first_nonzero_elem;
	int first_diag_index, first_in_row_index, first_in_column_index;
	
	//di[0] = di[0]  // it don't change.

	for (i = 1; i < n; i++) {	// loop through ia
		column = i - ia[i + 1] + ia[i];
		for (j = ia[i]; j < ia[i + 1]; j++, column++) {	// loop through al & au
			sumL = sumU = 0;
			amount_elems = std::min(j - ia[i], ia[column + 1] - ia[column]);
			first_diag_index = std::max(column - j + ia[i], column - ia[column + 1] + ia[column]);
			ia[column + 1] - ia[column] < j - ia[i] ? offset = j - ia[i] - ia[column + 1] + ia[column] : offset = 0;
		
			first_in_row_index = ia[i] + offset;					
			first_in_column_index = ia[column + 1] - amount_elems;

			for (k = 0; k < amount_elems; k++) {
				sumL += al[first_in_row_index]   * au[first_in_column_index]   * di[first_diag_index];
				sumU += au[first_in_row_index++] * al[first_in_column_index++] * di[first_diag_index++]; // ?
				//printf_s("sumL: %f %f %f\n", al[ia[i] + k + offset], di[first_diag_index + k], au[ia[column + 1] - amount_elems + k]);
				//printf_s("sumU: %f %f %f\n", al[ia[column + 1] - amount_elems + k], di[first_diag_index + k], au[ia[i] + k + offset]);
			}

			al[j] = (al[j] - sumL) / di[column];
			au[j] = (au[j] - sumU) / di[column];
		}
	
		// calculate diagonal elem:
		sumL = 0;

		for (k = ia[i], j = i - ia[i + 1] + ia[i]; k < ia[i + 1]; k++, j++) {
		//	printf_s("diag: %f %f %f\n", al[ia[i] + k], di[i - ia[i + 1] + ia[i] + k], au[ia[i] + k]);
			sumL += al[k] * di[j] * au[k]; 
		}

		if (abs((di[i] - sumL)) < abs(di[0])*eps) {
			printf("LDU decomposition don't possible, zero on diagonal");
			_getch();
			return;
		}
		else {
			di[i] -= sumL;
		}

	}
}


Full_matrix Matrix::in_full() {
	Full_matrix M(n);
	int column;

	for (int i = 0; i < n; i++) {
		column = i - ia[i + 1] + ia[i];
		for (int j = ia[i]; j < ia[i + 1]; j++, column++) {
			M.set_value(i, column, al[j]);
			M.set_value(column, i, au[j]);
		}
	}

	for (int i = 0; i < n; i++) {
		M.set_value(i, i, di[i]);
	}

	return M;
}

std::vector<mytype> Matrix::multipl_with_vector(std::vector<mytype>& x) {
	std::vector<mytype> b(n);
	int column;

	for (int i = 0; i < n; i++) {
		column = i - ia[i + 1] + ia[i];
		for (int j = ia[i]; j < ia[i + 1]; j++, column++) {
			b[i] += al[j] * x[column];
			b[column] += au[j] * x[i];
		}
	}

	for (int i = 0; i < n; i++) {
		b[i] += di[i] * x[i];
	}
	
	return b;
}

std::vector<mytype> Matrix::solve_system(std::vector<mytype>& b , std::vector<mytype>& y) {
	mytype2 sum;
	int column;

	for (int i = 0; i < n; i++) {			// Calculate Y
		sum = 0;
		column = i - ia[i + 1] + ia[i];
		for (int j = ia[i]; j < ia[i + 1]; j++, column++) {
			sum += y[column] * al[j];
		}
		y[i] = b[i] - sum;
	}

	std::vector<mytype> *z = &b;

	for (int i = 0; i < n; i++) {			// Calculate Z
		(*z)[i] = y[i] / di[i];
	}

	std::vector<mytype> *x = &y;

	for (int i = n - 1; i >= 0; i--) {		// Calculate X
		(*x)[i] = (*z)[i];
		column = i - ia[i + 1] + ia[i];
		for (int j = ia[i]; j < ia[i + 1]; j++, column++) {
			(*z)[column] -= (*x)[i] * au[j];
		}
	}

	return *x;
}

void Matrix::research1(std::vector<mytype>& x, std::vector<mytype>& b) {
	FILE *ans;
	FILE *err;
	fopen_s(&ans, "ans.txt", "w");
	fopen_s(&err, "err.txt", "w");

	for (int i = 0; i < 17; i++) {
		load();
		x = load_vector(n, "x.txt", x);
		di[0] += pow(10, -i);

		b = multipl_with_vector(x);
		LDU_decompose();
		x = solve_system(b, x);

		//Full_matrix Mf = in_full();
		//x = Mf.method_Gaussa(x, b);

		fprintf_s(ans, "i = %d: \n", i);
		fprintf_s(err, "i = %d: \n", i);

		for (int k = 0; k < n; k++) fprintf_s(ans, "%f\n", x[k]);
		for (int k = 0; k < n; k++) fprintf_s(err, "%e\n", (k + 1) - x[k]);

		fprintf_s(ans, "\n");
		fprintf_s(err, "\n");
	}

	fclose(ans);
	fclose(err);
}

void Matrix::change(int n1, int m1, std::vector<int>& ia1, std::vector<double>& al1, std::vector<double>& au1, std::vector<double>& di1) {
	n = n1; m = m1;
	//ia.resize(ia1.size());
	ia = ia1;

	al.resize(al1.size());
	for (int i = 0; i < al1.size(); i++) al[i] = (mytype)al1[i];

	au.resize(au1.size());
	for (int i = 0; i < au1.size(); i++) au[i] = (mytype)au1[i];

	di.resize(di1.size());
	for (int i = 0; i < di1.size(); i++) di[i] = (mytype)di1[i];
}

/*
void Matrix::research2(std::vector<mytype>& x, std::vector<mytype>& b) {
	FILE *ans;
	FILE *err;
	fopen_s(&ans, "ans.txt", "w");
	fopen_s(&err, "err.txt", "w");


	for (int i = 1; i < 17; i++) {		

		x.resize(i);
		for (int k = 0; k < i; k++)  x[k] = k + 1;

		//b = multipl_with_vector(x);
		//LDU_decompose();
		//x = solve_system(b, x);

		Full_matrix Mf = in_full();
		x = Mf.method_Gaussa(x, b);

		fprintf_s(ans, "i = %d: \n", i);
		fprintf_s(err, "i = %d: \n", i);

		for (int k = 0; k < n; k++) fprintf_s(ans, "%.15lf\n", x[k]);
		for (int k = 0; k < n; k++) fprintf_s(err, "%e\n", (k + 1) - x[k]);

		fprintf_s(ans, "\n");
		fprintf_s(err, "\n");
	}

	fclose(ans);
	fclose(err);
}
*/

/*
void Matrix::generate_Gilbert_matrix(int dim) {
	n = dim;

	ia.resize(n + 1);

	ia[0] = 0;
	for (int i = 1; i < n + 1; i++) ia[i] = ia[i - 1] + i - 1;

	m = ia[n];

	au.resize(m);
	al.resize(m);
	di.resize(n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < ia[i+1] - ia[i]; j++) {
			al[ia[i] + j] = au[ia[i] + j] = 1.0 / (i + j + 1);
		}
	}

	for (int i = 0; i < n; i++) di[i] = 1.0 / (i * 2 + 1);
}
*/

void Full_matrix::set_value(int row, int column, mytype value) {
	elems[row][column] = value;
}

Full_matrix::Full_matrix(int dim) {
	elems = std::vector<std::vector<mytype>>(dim);
	for (int i = 0; i < dim; i++) {
		elems[i] = std::vector<mytype>(dim);
	}
}

void Full_matrix::print() {
	for (int i = 0; i < elems.size(); i++) {
		for (int j = 0; j < elems[i].size(); j++) {
			printf_s(/*frmt.c_str()*/ "%6.4lf ", elems[i][j]);
		}
		printf_s("\n");
	}
}

void out_vector(std::vector<mytype>& x) {
	FILE *in;
	fopen_s(&in, "x.txt", "w");
	for (int i = 0; i < x.size(); i++) {
		fprintf_s(in, frmt.c_str(), x[i]);
	}
	fclose(in);
}

Full_matrix::Full_matrix() {};
Full_matrix::~Full_matrix() {
	for (int i = 0; i < elems.size(); i++) {
		elems[i].clear();
	}
	elems.clear();
}

std::vector<mytype> Full_matrix::method_Gaussa(std::vector<mytype>& x, std::vector<mytype>& b) {
	int i, j, k, index_main_row = 0;
	double max, sum, coef;

	for (i = 0; i < elems.size(); i++) {
		max = std::abs(elems[i][i]);
		index_main_row = i;

		for  (j = i + 1; j < elems.size(); j++) {
			if (std::abs(elems[j][i]) > std::abs(max)) {
				index_main_row = j;
				max = elems[j][i];
			}
		}

		// rows swap:
		if (i != index_main_row) {
			std::swap(elems[i], elems[index_main_row]);
			std::swap(b[i], b[index_main_row]);
		}

		// column zeroing:
		for (j = i + 1; j < elems.size(); j++) {
			coef = elems[j][i] / elems[i][i];

			for (k = i; k < elems[i].size(); k++) {
				elems[j][k] -= elems[i][k] * coef;
			}

			b[j] -= b[i] * coef;
		}
	}

	// back 

	for (i = x.size() - 1; i >= 0; i--) {
		sum = 0;

		for (j = x.size() - 1; j > i; j--) {
			sum += elems[i][j] * x[j];
		}

		x[i] = (b[i] - sum) / elems[i][i];
	}

	return x;
}