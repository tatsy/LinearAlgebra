#include <vector>
using namespace std;

#define __MAT64F_EXPORT__
#include "Matrix64f.h"
using namespace funopt;

#include "MTRand.h"
#include "Vector64f.h"
#include "funopt_macros.h"

Matrix64f::Matrix64f() :
    nrows(0),
    ncols(0),            
    data(0)
{
}

Matrix64f::Matrix64f(int rows, int cols) :
    nrows(rows),
    ncols(cols),
    data(0)
{
    data = new double[rows * cols];
}

Matrix64f::Matrix64f(double* data_, int rows, int cols) :
    nrows(rows),
    ncols(cols),
    data(0)
{
    data = new double[nrows * ncols];
    memcpy(data, data_, sizeof(double) * nrows * ncols);
}

Matrix64f::Matrix64f(const Matrix64f& m) :
    nrows(m.nrows),
    ncols(m.ncols),
    data(0)
{
    data = new double[nrows * ncols];
    memcpy(data, m.data, sizeof(double) * nrows * ncols);
}

Matrix64f::Matrix64f(const Vector64f& v) :
	nrows(v.ndim),
	ncols(1),
	data(0)
{
	data = new double[nrows * ncols];
	memcpy(data, v.data, sizeof(double) * nrows * ncols);
}

// 単位行列
Matrix64f Matrix64f::eye(int n) {
	Matrix64f I(n, n);
	memset(I.data, 0, sizeof(double) * n * n);
	for(int i=0; i<n; i++) {
		I(i, i) = 1.0;
	}
	return I;
}

// 零行列
Matrix64f Matrix64f::zeros(int rows, int cols) {
	Matrix64f Z(rows, cols);
	memset(Z.data, 0, sizeof(double) * rows * cols);
	return Z;
}

// 乱数行列
Matrix64f Matrix64f::rand(int rows, int cols)
{
	MTRand rng;
	Matrix64f ret(rows, cols);
	for(int i=0; i<rows; i++) {
		for(int j=0; j<cols; j++) {
			ret(i, j) = rng.randReal();
		}
	}
	return ret;
}

// 対角行列
Matrix64f Matrix64f::diag(vector<double>& v)
{
    const int n = (int)v.size();
    Matrix64f D(n, n);
    memset(D.data, 0, sizeof(double) * n * n);
    for(int i=0; i<n; i++) {
        D(i, i) = v[i];
    }
    return D;
}

// デストラクタ
Matrix64f::~Matrix64f()
{
    delete[] data;
}

Matrix64f& Matrix64f::operator=(const Matrix64f& m)
{
    delete[] data;

    nrows    = m.nrows;
    ncols    = m.ncols;
    data     = new double[nrows * ncols];
    memcpy(data, m.data, sizeof(double) * nrows * ncols);
    return *this;
}

inline double& Matrix64f::operator()(int i, int j) 
{
    return data[i * ncols + j];
}

inline double Matrix64f::operator()(int i, int j) const 
{
    return data[i * ncols + j];
}

inline int Matrix64f::rows() const 
{
    return nrows;
}

inline int Matrix64f::cols() const 
{
    return ncols;
}

Matrix64f Matrix64f::submat(int i, int j, int rows, int cols) const
{
    massert(i + rows <= nrows && j + cols <= ncols, "Too large sub-matrix size !!");

    Matrix64f A(rows, cols);
    for(int ii=0; ii<rows; ii++) {
        for(int jj=0; jj<cols; jj++) {
            A(ii, jj) = (*this)(i+ii, j+jj);
        }
    }
    return A;
}

Matrix64f Matrix64f::trans() const {
	Matrix64f ret(ncols, nrows);
	for(int i=0; i<ncols; i++) {
		for(int j=0; j<nrows; j++) {
			ret(i, j) = (*this)(j, i);
		}
	}
	return ret;
}

double Matrix64f::det() const {
	massert(nrows == ncols, "Matrix is not square. Cannot factorize.");
	Matrix64f LU;
	int* order = new int[nrows];
	factor_lu(LU, order);
	delete[] order;

	double ret = 1.0;
	for(int i=0; i<nrows; i++) {
		ret *= LU(i, i);
	}
	return ret;
}

Matrix64f Matrix64f::inv() const {
	massert(nrows == ncols, "Matrix is not square. Cannot invert.");

	int n = nrows;
	Matrix64f I = Matrix64f::eye(n);
	Matrix64f B;
	solve_lu(I, B); 
	return B;
}

Matrix64f Matrix64f::solve(const Matrix64f& b, SolverType solver_type)
{
    Matrix64f x;
    switch(solver_type) {
    case SOLVER_LU:
	    solve_lu(b, x);
        break;

    case SOLVER_QR:
		solve_qr(b, x);
        break;

    case SOLVER_CG:
        solve_cg(b, x);
        break;
	}
    return x;
}

ostream& operator<<(ostream& os, const Matrix64f& m)
{
    for(int i=0; i<m.rows(); i++) {
        os << (i == 0 ? "[" : " ");
        os << "[ ";
        for(int j=0; j<m.cols(); j++) {
            os << m(i, j) << (j == m.cols()-1 ? " ]" : ", ");
        }
        os << (i == m.rows()-1 ? "]" : "\n");
    }
    return os;
}
