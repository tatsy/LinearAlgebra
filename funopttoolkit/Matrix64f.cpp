#define __EXPORT__
#include "Matrix64f.h"
using namespace funopt;

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

// デストラクタ
Matrix64f::~Matrix64f()
{
    delete[] data;
}

Matrix64f& Matrix64f::operator=(const Matrix64f& m)
{
    delete[] data;

    nrows = m.nrows;
    ncols = m.ncols;
    data  = new double[nrows * ncols];
    memcpy(data, m.data, sizeof(double) * nrows * ncols);

    return *this;
}

double& Matrix64f::operator()(int i, int j) {
    return data[i * ncols + j];
}

double Matrix64f::operator()(int i, int j) const {
    return data[i * ncols + j];
}

Matrix64f Matrix64f::operator+(const Matrix64f& m) const {
	massert(nrows == m.nrows && ncols == m.ncols, "Matrix size is invalid");
	Matrix64f ret(nrows, ncols);
	int n = nrows * ncols;
	for(int i=0; i<n; i++) {
		ret.data[i] = data[i] + m.data[i];
	}
	return ret;
}

Matrix64f Matrix64f::operator-(const Matrix64f& m) const {
	massert(nrows == m.nrows && ncols == m.ncols, "Matrix size is invalid");
	Matrix64f ret(nrows, ncols);
	int n = nrows * ncols;
	for(int i=0; i<n; i++) {
		ret.data[i] = data[i] - m.data[i];
	}
	return ret;
}


Vector64f Matrix64f::operator*(const Vector64f& v) const {
    Vector64f ret(v.ndim);
    for(int i=0; i<nrows; i++) {
        double val = 0.0;
        for(int j=0; j<ncols; j++) {
            val += (*this)(i, j) * v(j);
        }
		ret(i) = val;
    }
    return ret;
}

Matrix64f Matrix64f::operator*(const Matrix64f& m) const {
	massert(ncols == m.nrows, "Matrix size is invalid");

	Matrix64f ret(nrows, m.ncols);
	for(int i=0; i<nrows; i++) {
		for(int j=0; j<m.ncols; j++) {
			double val = 0.0;
			for(int k=0; k<ncols; k++) {
				val += (*this)(i, k) * m(k, j);
			}
			ret(i, j) = val;
		}
	}
	return ret;
}

Matrix64f Matrix64f::operator*(const double d) const {
	Matrix64f ret(nrows, ncols);
	int n = nrows * ncols;
	for(int i=0; i<n; i++) {
		ret.data[i] = data[i] * d;
	}
	return ret;
}

Matrix64f Matrix64f::operator/(const double d) const {
	Matrix64f ret(nrows, ncols);
	int n = nrows * ncols;
	for(int i=0; i<n; i++) {
		ret.data[i] = data[i] / d;
	}
	return ret;
}


int Matrix64f::rows() const {
    return nrows;
}

int Matrix64f::cols() const {
    return ncols;
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

double Matrix64f::norm() const {
	return sqrt(norm2());
}

double Matrix64f::norm2() const {
	double ret = 0.0;
	int n = nrows * ncols;
	for(int i=0; i<n; i++) {
		ret += data[i] * data[i];
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

Matrix64f Matrix64f::solve(Matrix64f& b, int factor_type) {
    Matrix64f x;
	if(factor_type == FUNOPT_FACTOR_LU) {
	    solve_lu(b, x);
	}
	else if(factor_type == FUNOPT_FACTOR_QR) {
		solve_qr(b, x);
	}
    return x;
}

