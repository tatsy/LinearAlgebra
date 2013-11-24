#include "Powell.h"
#include "Matrix64f.h"

namespace funopt {
    namespace nonlin {
        Powell::Powell()
        {
        }

        void Powell::solve(const funcNd& func, const Vector64f& x0, Vector64f& x_opt, const int maxiter, const double tol)
        {
            minimize(func, x0, maxiter, tol);
            x_opt = xmin;
        }

        double Powell::minimize(const funcNd& func, const Vector64f& x0, const int maxiter, const double tol)
        {
            static const double TINY = 1.0e-20;
            const int n = x0.dim();
           
            //　初期方向は各軸方向に設定
            Matrix64f emat = Matrix64f::eye(n);

            Vector64f x_cur = xmin = x0;
            Vector64f e_cur(n);
            double fret = func(x0);
            for(int it=0; it<maxiter; it++) {
                // 現在の各方向に順に進む
                // その時の最大減少量を保存しておく (終了判定に使う)
                double f0 = fret;
                double max_dec = 0.0;
                int    max_i   = 0;
                for(int i=0; i<n; i++) {
                    x_cur = xmin;
                    for(int d=0; d<n; d++) {
                        e_cur(d) = emat(d, i);
                    }
                    double fsave = fret;
                    fret = linmin(func, x_cur, e_cur, maxiter, tol);
                    if(f0 - fret > max_dec) {
                        max_dec = fsave - fret;
                        max_i   = i;
                    }
                }

                // 終了判定
                if(2.0 * (f0 - fret) <= tol * (abs(f0) + abs(fret)) + TINY) {
                    break;
                }

                // 外分点の評価値を調べる
                Vector64f xnew = 2.0 * x_cur - xmin;                
                Vector64f enew = x_cur - xmin;
                double fnew = func(xnew);
                if(f0 > fnew) {
                    double fac1 = f0 - fret - max_dec;
                    double fac2 = f0 - fnew;
                    double t = 2.0 * (f0 - 2.0 * fret + fnew) * fac1 * fac1 - max_dec * fac2 * fac2;
                    if(t < 0.0) {
                        fret = linmin(func, xnew, enew, maxiter, tol);
                        for(int j=0; j<n; j++) {
                            emat(j, max_i) = emat(j, n-1);
                            emat(j, n-1)   = xmin(j);
                        }                        
                    }
                }
            }
            return fret;
        }
    }
}
