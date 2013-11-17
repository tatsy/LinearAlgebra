#define __EXPORT__
#include "Vector64f.h"
#include "Matrix64f.h"
using namespace funopt;

static double matrix_norm_frobenius(const Matrix64f& A)
{
    int r = A.rows();
    int c = A.cols();
    double ret = 0.0;
    for(int i=0; i<r; i++) {
        for(int j=0; j<c; j++) {
            ret += A(i, j) * A(i, j);
        }
    }
    ret = sqrt(ret);
    return ret;
}

static double matrix_norm_spectral(const Matrix64f& A)
{
    const double tol = 1.0e-20;
    const int    dim = A.cols();
    
    Vector64f u, v(dim);
    double vn, ret, old = 0.0;

    v(0) = 1.0;
    vn = v.norm();
    for(int it=0; it<100; it++) {
        u = A * v;
        ret = v.dot(u);
        if(abs(old - ret) < tol) break;

        old = ret;
        vn = u.norm();
        v = u / vn;
    }
    return ret;
}

double Matrix64f::norm(MatrixNormType type) const
{
    double ret = 0.0;
    switch(type) {
    case MAT_NORM_FROBENIUS:
        ret = matrix_norm_frobenius(*this);
        break;

    case MAT_NORM_SPECTRAL:
        ret = matrix_norm_spectral(*this);
        break;
    }
    return ret;
}
