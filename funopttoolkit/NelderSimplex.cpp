#include "NelderSimplex.h"

#include "Matrix64f.h"
#include "Vector64f.h"

namespace funopt {
    namespace nonlin {
        NelderSimplex::NelderSimplex()
        {
        }

        void NelderSimplex::solve(const funcNd& func, const Vector64f& x0, Vector64f& x_opt, const int maxiter, const double tol)
        {
            const int n = x0.dim();

            // 点の情報を入れた行列
            // 各列が点の座標に対応
            Vector64f* X = new Vector64f[n+1];
            X[0] = x0;
            for(int j=1; j<=n; j++) {
                X[j] = Vector64f(n);
                for(int i=0; i<n; i++) {
                    if(j == i+1) {
                        X[j](i) = X[0](i) + 1.0;
                    } else {
                        X[j](i) = X[0](i);
                    }
                }
            }

            // 点における値を評価
            double* fval = new double[n+1];
            for(int j=0; j<=n; j++) {
                fval[j] = func(X[j]);
            }

            for(int it=0; it<maxiter; it++) {
                // 最大値を持つ点を見つける
                int ihi  = 0;
                int ilo  = 0;
                int inhi = 0;
                for(int j=0; j<=n; j++) {
                    if(fval[ihi] < fval[j]) {
                        inhi = ihi;
                        ihi  = j;
                    } else if(fval[inhi] < fval[j]) {
                        inhi = j;
                    }

                    if(fval[j] < fval[ilo]) {
                        ilo = j;
                    }
                }

                // 終了判定
                if(abs(fval[ihi] - fval[ilo]) < tol) break;

                // 最大値を持つ点をアップデート
                update(func, X, fval, ihi, -1.0); 

                // 結果が良ければもっと進む
                if(fval[ihi] <= fval[ilo]) {
                    update(func, X, fval, ihi, 2.0);
                }
                // あまり下がらなかったら少し戻る
                else if(fval[ihi] >= fval[inhi]) {
                    double fhi = fval[ihi];
                    update(func, X, fval, ihi, 0.5);

                    // それでもあまり下がらなければ
                    // 最良点に他の点を近づける
                    if(fhi <= fval[ihi]) {
                        for(int j=0; j<=n; j++) {
                            if(j == ilo) continue;
                            X[j] = (X[j] + X[ilo]) * 0.5;
                            fval[j] = func(X[j]);
                        }
                    }
                }
            }
            x_opt = X[0];

            delete[] X;
            delete[] fval;
        }

        void NelderSimplex::update(const funcNd& func, Vector64f* xs, double* fs, int ihi, double t)
        {
            const int n = xs[0].dim();
            Vector64f x_new = Vector64f::zeros(n);
            for(int i=0; i<n+1; i++) {
                if(i == ihi) continue;
                x_new += xs[i];
            }
            x_new   = t * xs[ihi] + (1.0 - t) * x_new / (double)n;
            double f_new = func(x_new);
            if(f_new < fs[ihi]) {
                xs[ihi] = x_new;
                fs[ihi] = f_new;
            }
        }
    }
}
