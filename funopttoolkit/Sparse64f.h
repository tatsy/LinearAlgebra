#ifndef _SPARSE_64F_
#define _SPARSE_64F_

#include <vector>
using namespace std;

#include "dll_macros.h"

namespace funopt {
	class DLL_EXPORT Triplet;

	class Sparse64f {
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

		static Sparse64f eye(int n);
		static Sparse64f diags(vector<double> entries, int n);
	};
}

#endif
