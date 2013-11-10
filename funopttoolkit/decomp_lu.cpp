#include <cstdio>
#include <cstdlib>

#include "funopt_macros.h"
#include "Matrix64f.h"
#include "Vector64f.h"

#include "linsolve.h"

void funopt::factor_lu(Matrix64f& A, Matrix64f& LU) {
    massert(A.rows() == A.cols(), "Matrix is not square. Cannot factorize.");
    
    int n = A.cols();
    LU = A;        
    for(int j=0; j<n; j++) {
        for(int i=1; i<=j; i++) {
            for(int k=0; k<i; k++) {
                LU(i, j) -= LU(i, k) * LU(k, j);
            }
        }

        double ujj = LU(j, j);
        for(int i=j+1; i<n; i++) {
            for(int k=0; k<j; k++) {
                LU(i, j) -= LU(i, k) * LU(k, j);
            }
            LU(i, j) /= ujj;
        }
    }
}

void funopt::solve_lu(Matrix64f& A, Vector64f& b, Vector64f& x) {
    Matrix64f LU;
    factor_lu(A, LU);

    int n = A.rows();
    x = b;
    // solve for L
    for(int j=0; j<n; j++) {
        for(int i=j+1; i<n; i++) {
            x(i) = x(i) - x(j) * LU(i, j);
        }
    }

    // solve for U
    for(int j=n-1; j>=0; j--) {
        x(j) = x(j) / LU(j, j);
        for(int i=0; i<j; i++) {
            x(i) = x(i) - x(j) * LU(i, j);
        }
    }
}
