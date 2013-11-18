#include <map>
#include <algorithm>
using namespace std;

#include "funopt_macros.h"

#define __EXPORT__
#include "Sparse64f.h"
using namespace funopt;

Sparse64f::Sparse64f() :
	nrows(0),
	ncols(0),
	data(0),
	row_idx(0),
	col_heads(0)
{
}

Sparse64f::Sparse64f(vector<Triplet> triplets, int rows, int cols) :
	nrows(rows),
	ncols(cols),
	nnz((int)triplets.size()),
	data(0),
	row_idx(0),
	col_heads(0)
{
	data = new double[nnz];
	row_idx = new int[nnz];
	col_heads = new int[ncols+1];

	sort(triplets.begin(), triplets.end());

	for(int i=0; i<nnz; i++) {
		data[i]    = triplets[i].val;
		row_idx[i] = triplets[i].i;
		if(i == 0) col_heads[0] = 0;
		else if(triplets[i].j != triplets[i-1].j) {
			col_heads[triplets[i].j] = i;
		}
	}
	col_heads[ncols] = nnz;
}

Sparse64f::Sparse64f(const Sparse64f& sp) :
	nrows(sp.nrows),
	ncols(sp.ncols),
	nnz(sp.nnz),
	data(0),
	row_idx(0),
	col_heads(0)
{
	data = new double[nnz];
	row_idx = new int[nnz];
	col_heads = new int[ncols+1];

	memcpy(data, sp.data, sizeof(double) * nnz);
	memcpy(row_idx, sp.row_idx, sizeof(int) * nnz);
	memcpy(col_heads, sp.col_heads, sizeof(int) * (ncols+1));
}

Sparse64f::~Sparse64f()
{
	delete[] data;
	delete[] row_idx;
	delete[] col_heads;
}

Sparse64f& Sparse64f::operator=(const Sparse64f& sp)
{
	delete[] data;
	delete[] row_idx;
	delete[] col_heads;

	nrows = sp.nrows;
	ncols = sp.ncols;
	nnz   = sp.nnz;

	data = new double[nnz];
	row_idx = new int[nnz];
	col_heads = new int[ncols+1];

	memcpy(data, sp.data, sizeof(double) * nnz);
	memcpy(row_idx, sp.row_idx, sizeof(int) * nnz);
	memcpy(col_heads, sp.col_heads, sizeof(int) * (ncols+1));

	return *this;
}

Sparse64f Sparse64f::operator*(const Sparse64f& sp) const
{
	massert(ncols == sp.nrows, "Matrix size is invalid");
	int rows = nrows;
	int cols = sp.ncols;
	int len  = ncols;

	map<int,double>::iterator it;
	map<int,double> col_elems;
	vector<Triplet> triplets;
	for(int k=0; k<len; k++) {
		col_elems.clear();
		int s1 = sp.col_heads[k];
		int t1 = sp.col_heads[k+1];
		for(int i=s1; i<t1; i++) {
			int r1  = row_idx[i];
			int s2 = col_heads[r1];
			int t2 = col_heads[r1+1];
			for(int j=s2; j<t2; j++) {
				int r2 = row_idx[j];
				col_elems[r2] += data[j] * sp.data[i];
			}
		}

		for(it=col_elems.begin(); it!=col_elems.end(); ++it) {
			triplets.push_back(Triplet(it->first, k, it->second));
		}
	}
	return Sparse64f(triplets, rows, cols);
}

double Sparse64f::operator()(int i, int j) const
{
	for(int k=col_heads[j]; k<col_heads[j+1]; k++) {
		if(row_idx[k] == i) {
			return data[k];
		}
	}
	return 0.0;
}

Sparse64f Sparse64f::eye(int n)
{
	vector<Triplet> triplets(n);
	for(int i=0; i<n; i++) {
		triplets[i] = Triplet(i, i, 1.0);
	}
	return Sparse64f(triplets, n, n);
}

Sparse64f Sparse64f::diags(vector<double> entries, int n)
{
	vector<Triplet> triplets(n);
	for(int i=0; i<n; i++) {
		triplets[i] = Triplet(i, i, entries[i]);
	}
	return Sparse64f(triplets, n, n);
}
