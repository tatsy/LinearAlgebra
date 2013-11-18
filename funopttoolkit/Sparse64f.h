#ifndef _SPARSE_64F_
#define _SPARSE_64F_

#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

#include "dll_macros.h"
#include "Triplet.h"

namespace funopt {
	class Triplet;

	class DLL_EXPORT Sparse64f {
	private:
		int nrows, ncols, nnz;
		double* data;
		int*    row_idx;
		int*    col_heads;

	public:
		Sparse64f();
		Sparse64f(vector<Triplet> triplets, int rows, int cols);
		Sparse64f(const Sparse64f& sp);
		~Sparse64f();

		Sparse64f& operator=(const Sparse64f& sp);
		Sparse64f  operator*(const Sparse64f& sp) const;
		double     operator()(int i, int j) const;

		static Sparse64f eye(int n);
		static Sparse64f diags(vector<double> entries, int n);
	};
}

#endif
