#include <vector>
#include <algorithm>
using namespace std;

#include "Vector64f.h"

#define __MAT64F_EXPORT__
#include "funopt_macros.h"
#include "Matrix64f.h"
using namespace funopt;

class abs_comparator {
public:
    abs_comparator() {}
    bool operator()(const double& a, const double& b) const 
    {
        return abs(a) > abs(b);
    }
};

// 逆べき乗法による固有値の計算
static void eig_invpower(const Matrix64f& A, Matrix64f& D, Matrix64f& U) 
{
    massert(A.rows() == A.cols(), "matrix is not square");
    const int n = A.rows();

    Matrix64f M = A;
        D = Matrix64f::zeros(n, n);
    U = Matrix64f::zeros(n, n);

    for(int i=0; i<n; i++) {
        Matrix64f b = Matrix64f::rand(n, 1);
        double nvold = b.norm();
        double s = 0.0;
        for(int it=0; it<1000; it++) {
            s = ((b.trans() * M * b) / (b.trans() * b)(0, 0))(0, 0);
            Matrix64f Mt = M - Matrix64f::eye(n) * s;
            Matrix64f bt = Mt.solve(b);

            double nv = bt.norm();
            b = bt / nv;
            if(abs(nv - nvold) < 1.0e-20) {
                break;
            }
            nvold = nv;
            cout << s << endl;
        }
        cout << s << endl;
        cout << b << endl;

        for(int j=0; j<n; j++) {
            U(j, i) = b(j, 0);
        }
        M = M - b * b.trans() * s;
        D(i, i) = s;
    }
}

// Power法による固有値の計算
static void eig_power(const Matrix64f& A, Matrix64f& D, Matrix64f& U) 
{
    massert(A.rows() == A.cols(), "matrix is not square");
    const int n = A.rows();

    Matrix64f M = A;
    D = Matrix64f::zeros(n, n);
    U = Matrix64f::zeros(n, n);
    vector<Matrix64f> us;

    for(int i=0; i<n; i++) {
        Matrix64f b = Matrix64f::rand(n, 1);
        double nvold = b.norm();
        for(int it=0; it<1000; it++) {
            Matrix64f bt = M * b;
            for(int k=0; k<us.size(); k++) {
                bt = bt - D(k, k) * (us[k] * (us[k].trans() * b));
            }

            double nv = bt.norm();
            b = bt / nv;
            if(abs(nv - nvold) < 1.0e-20) {
                break;
            }
            nvold = nv;
        }

        D(i, i) = (b.trans() * M * b)(0, 0);
        for(int j=0; j<n; j++) {
            U(j, i) = b(j, 0);
        }
        us.push_back(b);
        //M = M - b * b.trans() * D(i, i);
    }
}

// ギブンス回転を用いて三重対角行列を対角化
static void givens_rotate(const Matrix64f& A, Matrix64f& Q, Matrix64f& R)
{
    massert(A.rows() == A.cols(), "matrix is not square");
    const int n = A.rows();
    Matrix64f I = Matrix64f::eye(n);

    R = A;
    Q = I;
    for(int k=0; k<n-1; k++) {
        Matrix64f G = I;
        double app = R(k, k);
        double aqp = R(k+1, k);
        double r   = sqrt(app * app + aqp * aqp);
        double c   =  app / r;
        double s   = -aqp / r;
        G(k,   k)   = c;
        G(k,   k+1) = -s;
        G(k+1, k)   = s;
        G(k+1, k+1) = c;
        R = G * R;
        swap(G(k, k+1), G(k+1, k));
        Q = Q * G;
    }
}

// 行列Aを三重対角行列Tと相似変換Vに分解
static void tridiag(const Matrix64f& A, Matrix64f& T, Matrix64f& V)
{
	massert(A.rows() == A.cols(), "matrix is not square");
	int n = A.rows();
    T = A;
    V = Matrix64f::eye(n);
	double nx = 0.0;
	for(int k=0; k<n-2; k++) {
        Matrix64f xt = Matrix64f::zeros(n, 1);
        for(int i=k+1; i<n; i++) {
            xt(i, 0) = T(i, k);
        }

        double xi = 0.0;
        double nx = 0.0;
        for(int i=k+2; i<n; i++) {
            xi += xt(i, 0) * xt(i, 0);
        }
        nx = sqrt(xi + xt(k+1, 0) * xt(k+1, 0));
        xi = xi / (abs(xt(k+1, 0)) + nx);


        Matrix64f u = xt;
        double sgn = xt(k+1, 0) > 0.0 ? 1.0 : -1.0;
        u(k+1,0) = -sgn * xi;

        double nu = 0.0;
        for(int i=k+1; i<n; i++) {
            nu += u(i, 0) * u(i, 0);
        }

        Matrix64f P = Matrix64f::eye(n) - 2.0 * u * u.trans() / nu;
        T = P * T * P;
        V = V * P;
	}
}

// 固有値分解
void Matrix64f::eig(Matrix64f& val, Matrix64f& vec) const
{
	massert(nrows == ncols, "Matrix is not square.");
    eig_power(*this, val, vec);
    cout << val << endl;

    double tol = 1.0e-20;
    int n = nrows;

    Matrix64f I;
    I = Matrix64f::eye(n);

    Matrix64f T, V;
	tridiag(*this, T, V);

    // シフト付きQR法
    Matrix64f M = T;
    val = Matrix64f::zeros(n, n);
    vec = I;
    Matrix64f Q, R;
    for(int it=0; it<1000; it++) {
        double a11 = M(n-2,n-2);
        double a12 = M(n-2,n-1);
        double a21 = M(n-1,n-2);
        double a22 = M(n-1,n-1);
        double a = a11 + a22;
        double b = a11 * a22 - a12 * a21;
        double x1 = (a + sqrt(a * a - 4.0 * b)) / 2.0;
        double x2 = (a - sqrt(a * a - 4.0 * b)) / 2.0;
        double mu = abs(x1 - a22) < abs(x2 - a22) ? x1 : x2;
        M = M - I * mu;
        givens_rotate(M, Q, R);
        M = R * Q + I * mu;

        if(abs(M(n-1, n-2)) < tol) {
            val(n-1, n-1) = M(n-1, n-1);
            M = M.submat(0, 0, n-1, n-1);
            n = n - 1;
            I = Matrix64f::eye(n);
            if(n == 1) {
                val(0, 0) = M(0, 0);
                break;
            }
        }
    }
    n = val.rows();

    vector<double> v(n);
    for(int i=0; i<n; i++) v[i] = val(i, i);
    sort(v.begin(), v.end(), abs_comparator());
    val = Matrix64f::diag(v);

    // 三重対角行列から固有ベクトルを求める
    vec = Matrix64f(n, n);
    for(int j=0; j<n; j++) {
        double l = val(j, j);

        vec(0, j) = 1.0;
        vec(1, j) = -(T(0, 0) - l) / T(0, 1);
        for(int k=2; k<n; k++) {
            vec(k, j) = - (vec(k-2, j) * T(k-1, k-2) + vec(k-1, j) * (T(k-1, k-1) - l)) / T(k-1, k); 
        }
    }
    vec = V * vec;

    vec.factor_qr(V, R);
    vec = V;
    /*
    for(int j=0; j<n; j++) {
        double nv = 0.0;
        for(int i=0; i<n; i++) {
            nv += vec(i, j) * vec(i, j);
        }
        nv = sqrt(nv);
        for(int i=0; i<n; i++) {
            vec(i, j) /= nv;
        }
    }*/
}


