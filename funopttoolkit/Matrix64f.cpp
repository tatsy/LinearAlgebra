#define __EXPORT__
#include "Matrix64f.h"
using namespace funopt;

#include "Vector64f.h"
#include "linsolve.h"

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

// デストラクタ
Matrix64f::~Matrix64f()
{
    delete[] data;
}

Matrix64f& Matrix64f::operator=(const Matrix64f& m)
{
    if(data != 0) {
        delete[] data;
    }

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

Vector64f Matrix64f::operator*(const Vector64f& v) {
    Vector64f ret(v.ndim);
    for(int i=0; i<nrows; i++) {
        ret(i) = 0.0;
        for(int j=0; j<ncols; j++) {
            ret(i) += (*this)(i, j) * v(j);
        }
    }
    return ret;
}

int Matrix64f::rows() const {
    return nrows;
}

int Matrix64f::cols() const {
    return ncols;
}

Vector64f Matrix64f::solve(Vector64f& b, int decomp_type) {
    Vector64f x;
    solve_lu((*this), b, x);
    return x;
}

// 固有値を求める
void Matrix64f::eig(Matrix64f& eval, Matrix64f& evec) const {
}

