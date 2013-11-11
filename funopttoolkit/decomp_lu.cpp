#include <cstdio>
#include <cstdlib>

#define __EXPORT__
#include "funopt_macros.h"
#include "Matrix64f.h"
#include "Vector64f.h"

using namespace funopt;

void Matrix64f::factor_lu(Matrix64f& LU, int* order) const {
    massert(nrows == ncols, "Matrix is not square. Cannot factorize.");
    
    int n = ncols;
    LU = (*this);
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
            for(int j=0; j<n; j++) {
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
}

void Matrix64f::solve_lu(Matrix64f& b, Matrix64f& x) const {
	massert(nrows == ncols, "Matrix is not square. Cannot factorize.");
	massert(ncols == b.nrows, "Matrix size is invalid");

    int m = nrows;
	int n = b.ncols;
    
	// LU分解
	Matrix64f LU;
    int* order = new int[m];
    factor_lu(LU, order);

	// 行列式の計算
	double d = 1.0;
	for(int i=0; i<m; i++) {
		d *= LU(i, i);
	}
	massert(d != 0.0, "Matrix is singular. Cannot solve.");

	// ピボットに従って要素を入れ替える
    x = Matrix64f(b.nrows, b.ncols);
    for(int i=0; i<m; i++) {
		for(int j=0; j<n; j++) {
	        x(i, j) = b(order[i], j);
		}
    }

    // solve for L
    for(int j=0; j<m; j++) {
        for(int i=j+1; i<m; i++) {
			for(int k=0; k<n; k++) {
	            x(i, k) = x(i, k) - x(j, k) * LU(i, j);
			}
        }
    }

    // solve for U
    for(int j=m-1; j>=0; j--) {
		for(int k=0; k<n; k++) {
			x(j, k) = x(j, k) / LU(j, j);
			for(int i=0; i<j; i++) {
				x(i, k) = x(i, k) - x(j, k) * LU(i, j);
			}
		}
    }

    delete[] order;
}
