#include "Vector64f.h"
#include "Matrix64f.h"
using namespace funopt;

double matrix_norm_frobenius(const Matrix64f& A)
{
    int r = A.rows();
    int c = A.cols();
    double ret = 0.0;
    for(int i=0; i<r; i++) {
        for(int j=0; j<c; j++) {
            ret += A(i, j);
        }
    }
    ret = sqrt(ret);
    return ret;
}

double matrix_norm_spectral(const Matrix64f& A)
{
    int dim = A.cols();
    Vector64f v(dim);
    v(0) = 1.0;
    for(int it=0; it<100; it++) {
        v = A * v;
        v = v / v.norm();
    }

    double ret = v.dot(A * v) / v.norm();
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
