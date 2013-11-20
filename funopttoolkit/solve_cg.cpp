#include "funopt_macros.h"

#define __MAT64F_EXPORT__
#include "Matrix64f.h"
#include "Vector64f.h"
using namespace funopt;

static void solve_cg_vec (const Matrix64f& A, const Vector64f& b, Vector64f& x)
{
    const int nrows = A.rows();
    const int ncols = A.cols();
    x = Vector64f::rand(ncols);
    Vector64f c = b - A * x;
    Vector64f r = A.trans() * c;
    Vector64f p = c;
    double c2 = c.dot(c);

    for(int k=0; k<nrows; k++) {
        double c1 = c2;
        Vector64f q = A * p;
        Vector64f s = A.trans() * q;
        double alph = p.dot(r) / q.dot(q);
        x = x + alph * p;
        c = c - alph * q;
        r = r - alph * s;
        c2 = c.dot(c);
        double beta = -r.dot(s) / q.dot(q);
        p = r + beta * p;
        if(c1 <= c2) {
            break;
        }
    }
}

void Matrix64f::solve_cg(const Matrix64f& b, Matrix64f& x) const
{
    massert(nrows == ncols && nrows == b.nrows, "Matrix and vector size is invalid");

    // ‘Oˆ—
    Matrix64f A = *this;
    double* d = new double[nrows];
    for(int i=0; i<nrows; i++) {
        d[i] = 1.0 / sqrt(A(i, i));
    }

    for(int i=0; i<nrows; i++) {
        for(int j=0; j<ncols; j++) {
            A(i, j) *= d[i] * d[j];
        }        
    }

    // CG–@‚Å‰ð‚­
    x = Matrix64f(ncols, b.ncols);
    for(int i=0; i<b.ncols; i++) {
        Vector64f v(b.nrows);
        for(int j=0; j<b.nrows; j++) {
            v(j) = b(j, i) * d[j];
        }
        Vector64f u;
        solve_cg_vec(A, v, u);
        for(int j=0; j<b.nrows; j++) {
            x(j, i) = u(j) * d[j];
        }
    }

    delete[] d;
}

