#define __MAT64F_EXPORT__
#include "funopt_macros.h"
#include "Matrix64f.h"
using namespace funopt;

static void tridiag(Matrix64f& A)
{
	massert(A.rows() == A.cols(), "matrix is not square");

	int n = A.rows();
	double nx = 0.0;
	for(int k=0; k<n-1; k++) {
		double a = 0.0;
		for(int i=k+1; i<n; i++) {
			a += A(i, k) * A(i, k);
		}

		double sgn = A(k+1, k) > 0.0 ? 1.0 : -1.0;
		a = -sgn * sqrt(a);
		double r = sqrt(0.5 * a * (a - A(k+1, k)));

		Matrix64f v(n-k, 1);
		v(0, 0) = 0.0;
		v(1, 0) = (A(k+1, k) - a) / (2.0 * r);

		for(int i=2; i<n-k; i++) {
			v(i, 0) = A(k+i, k) / (2.0 * r);
		}

		Matrix64f subP = Matrix64f::eye(n - k) - (v * v.trans()) * 2.0;
		Matrix64f P = Matrix64f::eye(n);
		for(int i=0; i<n-k; i++) {
			for(int j=0; j<n-k; j++) {
				P(k+i, k+j) = subP(i, j);
			}
		}

		A = P * A * P;
	}
}

void Matrix64f::eig(Matrix64f& val, Matrix64f& vec) const
{
	massert(nrows == ncols, "Matrix is not square.");
    double tol = 1.0e-12;
    int n = nrows;

    Matrix64f Q, R, I;
    I = Matrix64f::eye(n);

    val = *this;
	tridiag(val);
	cout << val << endl << endl;

    vec = I;
    for(int i=0; i<1000; i++) {
        double a11 = val(n-2,n-2);
        double a12 = val(n-2,n-1);
        double a21 = val(n-1,n-2);
        double a22 = val(n-1,n-1);
        double a = a11 + a22;
        double b = a11 * a22 - a12 * a21;
        double x1 = (a + sqrt(a * a - 4.0 * b)) / 2.0;
        double x2 = (a - sqrt(a * a - 4.0 * b)) / 2.0;
        double mu = abs(x1 - a22) < abs(x2 - a22) ? x1 : x2;
        val = val - I * mu;
        val.factor_qr(Q, R);
        val = R * Q + I * mu;
        vec = vec * Q;

        bool is_end = true;
        for(int i=0; i<n-1; i++) {
            if(abs(val(i, i+1)) > tol || abs(val(i, i+1)) > tol) {
                is_end = false;
                break;
            }
        }

        if(is_end) break;
    }
}
