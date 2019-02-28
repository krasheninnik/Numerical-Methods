#pragma once
#include "vector"
#include "string"


//typedef float mytype;
typedef double mytype;
typedef double mytype2;

class Full_matrix{
public:
	Full_matrix();
	Full_matrix(int size);
	~Full_matrix();
	void set_value(int row, int column, mytype value);
	void print();
	std::vector<mytype>  method_Gaussa(std::vector<mytype>&, std::vector<mytype>&);

private:
	std::vector<std::vector<mytype>> elems;
};

// class of profile format matrix
class Matrix {
public:
	Matrix();
	~Matrix();

	Full_matrix in_full();
	int get_dim();
	void load();
	void LDU_decompose();
	std::vector<mytype> multipl_with_vector(std::vector<mytype>&);
	std::vector<mytype> solve_system(std::vector<mytype>&, std::vector<mytype>&);

	//void calcDiag(int);
	void research1(std::vector<mytype>&, std::vector<mytype>&);
	void research2(std::vector<mytype>&, std::vector<mytype>&);
	void change(int, int, std::vector<int>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);


private:
	int n;						// dimension of matrix
	int m;						// number of elems in up/low trinagle

	std::vector<int> ia;		// info about matrix
	std::vector<mytype> au;		// elem of upper triangle
	std::vector<mytype> al;		// elem of low triangle
	std::vector<mytype> di;		// diagonal elem
};
