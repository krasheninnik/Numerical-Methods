#include "matrix.h"

std::vector<mytype> load_vector(int, std::string, std::vector<mytype>&);
void out_vector(std::vector<mytype>&);

void main() {
	Matrix M;
	Full_matrix Mf;

	M.load();
	std::vector<mytype> x(M.get_dim()), b(M.get_dim());

	//b = load_vector(M.get_dim(), "b.txt", b);
	x = load_vector(M.get_dim(), "x.txt", x);
	b = M.multipl_with_vector(x);

	M.LDU_decompose();	
	x = M.solve_system(b, x);
	out_vector(x);
}