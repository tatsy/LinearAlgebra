#include "funopt_macros.h"

#define __MAT64F_EXPORT__
#include "Vector64f.h"
#include "Matrix64f.h"
using namespace funopt;

Matrix64f& Matrix64f::operator+=(const Matrix64f& A)
{
    massert(nrows == A.nrows && ncols == A.ncols, "Matrix size is different !!");
    int n = nrows * ncols;
    for(int i=0; i<n; i++) {
        data[i] += A.data[i];
    }
    return *this;
}

Matrix64f& Matrix64f::operator-=(const Matrix64f& A)
{
    massert(nrows == A.nrows && ncols == A.ncols, "Matrix size is different !!");
    int n = nrows * ncols;
    for(int i=0; i<n; i++) {
        data[i] -= A.data[i];
    }
    return *this;
}

Matrix64f operator+(const Matrix64f& A, const Matrix64f& B) 
{
	massert(A.rows() == B.rows() && A.cols() == B.cols(), "Matrix size is invalid");
    Matrix64f C = A;
    C += B;
    return C;
}

Matrix64f operator-(const Matrix64f& A, const Matrix64f& B) 
{
	massert(A.rows() == B.rows() && A.cols() == B.cols(), "Matrix size is invalid");
    Matrix64f C = A;
    C -= B;
    return C;
}

Vector64f operator*(const Matrix64f& A, const Vector64f& v)
{
    massert(A.cols() == v.dim(), "Matrix and vector size is invalid");
    int nrows = A.rows();
    int ncols = A.cols();
    Vector64f ret(nrows);
    for(int i=0; i<nrows; i++) {
        double val = 0.0;
        for(int j=0; j<ncols; j++) {
            val += A(i, j) * v(j);
        }
		ret(i) = val;
    }
    return ret;
}

Matrix64f operator*(const Matrix64f& A, const Matrix64f& B)
{
	massert(A.cols() == B.rows(), "Matrix size is invalid");
    const int m     = A.cols();
    const int nrows = A.rows();
    const int ncols = B.cols();
    Matrix64f C(nrows, ncols);
    for(int i=0; i<nrows; i++) {
        for(int j=0; j<ncols; j++) {
            double val = 0.0;
            for(int k=0; k<m; k++) {
                val += A(i, k) * B(k, j);
            }
            C(i, j) = val;
        }
    }
    return C;
}

Matrix64f operator*(const Matrix64f& A, const double d)
{
    const int nrows = A.rows();
    const int ncols = A.cols();
	Matrix64f C = A;
	for(int i=0; i<nrows; i++) {
        for(int j=0; j<ncols; j++) {
            C(i, j) *= d;
        }
	}
	return C;
}

Matrix64f operator*(double d, const Matrix64f& A)
{
    return A * d;
}

Matrix64f operator/(const Matrix64f& A, const double d)
{
    const int nrows = A.rows();
    const int ncols = A.cols();
	Matrix64f C = A;
	for(int i=0; i<nrows; i++) {
        for(int j=0; j<ncols; j++) {
            C(i, j) /= d;
        }
	}
	return C;
}

