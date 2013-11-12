#include <cstdio>
#include <cstdlib>

#define __EXPORT__
#include "funopt_macros.h"
#include "Matrix64f.h"
#include "Vector64f.h"
using namespace funopt;

void Matrix64f::factor_qr(Matrix64f& Q, Matrix64f& R) const {
	massert(nrows == ncols, "Matrix is not square. Cannot factorize.");

	int n = nrows;
	R = (*this);
	Q = Matrix64f::eye(n);
	for(int k=0; k<n-1; k++) {
		Matrix64f x(n - k, 1);
		for(int i=0; i<n - k; i++) {
			x(i, 0) = R(k+i, k);
		}

		x(0, 0) -= x.norm();
		Matrix64f subP = Matrix64f::eye(n - k) - (x * x.trans()) / (0.5 * x.norm2());
		Matrix64f P = Matrix64f::eye(n);
		for(int i=0; i<n-k; i++) {
			for(int j=0; j<n-k; j++) {
				P(k+i, k+j) = subP(i, j);
			}
		}

		R = P * R;
		Q = Q * P;
	}
}

void Matrix64f::solve_qr(Matrix64f& b, Matrix64f& x) const {
	massert(nrows == ncols, "Matrix is not square. Cannot factorize.");
	massert(ncols == b.nrows, "Matrix size is invalid");

    int m = nrows;
	int n = b.ncols;

	Matrix64f Q, R;
	factor_qr(Q, R);

	double d = 1.0;
	for(int i=0; i<n; i++) {
		d *= R(i, i);
	}
	massert(d != 0.0, "Matrix is sigular");

	x = b;
	x = Q.trans() * x;

	for(int j=m-1; j>=0; j--) {
		for(int k=0; k<n; k++) {
			x(j, k) = x(j, k) / R(j, j);
			for(int i=0; i<j; i++) {
				x(i, k) = x(i, k) - x(j, k) * R(i, j);
			}
		}
    }
}
