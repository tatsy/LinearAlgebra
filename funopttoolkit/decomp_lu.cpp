#include <cstdio>
#include <cstdlib>

#include "funopt_macros.h"
#include "Matrix64f.h"
#include "Vector64f.h"

#include "linsolve.h"

void funopt::factor_lu(Matrix64f& A, Matrix64f& LU, int* order) {
    massert(A.rows() == A.cols(), "Matrix is not square. Cannot factorize.");
    
    int n = A.cols();
    LU = A;
    for(int i=0; i<n; i++) order[i] = i;

    for(int k=0; k<n; k++) {
        // ピボット選択
        double maxval = 0.0;
        int    pivot  = k;
        for(int i=k; i<n; i++) {
            if(maxval < abs(LU(i, k))) {
                maxval = abs(LU(i, k));
                pivot  = i;
            }
        }

        // 行の入れ替え
        swap(order[k], order[pivot]);
        if(pivot != k) {
            for(int j=k; j<n; j++) {
                swap(LU(k, j), LU(pivot, j));
            }
        }

        // 要素の消去
        double iukk = 1.0 / LU(k, k);
        for(int i=k+1; i<n; i++) {
            LU(i, k) *= iukk;
            for(int j=k+1; j<n; j++) {
                LU(i, j) -= LU(i, k) * LU(k, j);
            }
        }
    }

    /*
    for(int j=0; j<n; j++) {
        for(int i=1; i<=j+1; i++) {
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
    */

    cout << LU << endl;
}

void funopt::solve_lu(Matrix64f& A, Vector64f& b, Vector64f& x) {
    
    int n = A.rows();
    Matrix64f LU;
    int* order = new int[n];
    factor_lu(A, LU, order);

    x = Vector64f(n);
    for(int i=0; i<n; i++) {
        printf("order[%d] = %d\n", i, order[i]);
        x(order[i]) = b(i);
    }
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

    delete[] order;
}
