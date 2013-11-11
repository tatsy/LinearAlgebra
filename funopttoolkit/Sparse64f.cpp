#include <algorithm>
using namespace std;

#include "Triplet.h"

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
	col_heads = new int[nrows];

	sort(triplets.begin(), triplets.end());

	for(int i=0; i<nnz; i++) {
		data[i]    = triplets[i].val;
		row_idx[i] = triplets[i].i;
		if(i == 0) col_heads[0] = 0;
		else if(triplets[i].j != triplets[i-1].j) {
			col_heads[triplets[i-1].j] = i;
		}
	}
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
	col_heads = new int[nrows];

	memcpy(data, sp.data, sizeof(double) * nnz);
	memcpy(row_idx, sp.row_idx, sizeof(int) * nnz);
	memcpy(col_heads, sp.col_heads, sizeof(int) * nrows);
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
	col_heads = new int[nrows];

	memcpy(data, sp.data, sizeof(double) * nnz);
	memcpy(row_idx, sp.row_idx, sizeof(int) * nnz);
	memcpy(col_heads, sp.col_heads, sizeof(int) * nrows);

	return *this;
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
