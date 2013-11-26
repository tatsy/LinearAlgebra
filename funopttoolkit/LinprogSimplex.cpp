#include "Matrix64f.h"
#include "Vector64f.h"
#include "LinprogSimplex.h"

namespace funopt {
    namespace linprog {
        Simplex::Simplex() :
            SolverBase()
        {
        }

        void Simplex::solve(const Vector64f& c_, const Matrix64f& A_, const Vector64f& b_)
        {
            Matrix64f A = A_;
            Vector64f c = c_;
            Vector64f b = b_;

            int  n   = A.rows();
            int  m   = A.cols();
            int* idx = new int[m];
            for(int i=0; i<m; i++) {
                idx[i] = i;
            }

            Vector64f x(n);
            for(;;) {
                // 基底に対応している行列成分を抜き出して逆行列を計算
                Matrix64f Ab    = A.submat(0, m-n, n, n);
                Matrix64f Abinv = Ab.inv();
                x = Abinv * b;
                
                Vector64f cb(n);
                for(int i=0; i<n; i++) cb(i) = c(m-n+i);
                Vector64f y = Abinv.trans() * cb;

                // 新しく入れる基底を選ぶ
                int    min_k = 0;
                double min_u = 1.0e20;
                for(int k=0; k<m-n; k++) {
                    Vector64f ak = A.col(k);
                    double u = c(k) - ak.dot(y);
                    if(min_u > u) {
                        min_u = u;
                        min_k = k;
                    }
                }
                if(min_u > 0.0) break;

                // 除外する基底を選ぶ
                Vector64f w = Abinv * A.col(min_k);
                int    min_i = 0;
                double min_r = 1.0e20;
                for(int i=0; i<n; i++) {
                    if(w(i) <= 0.0) continue;

                    double ratio = x(i) / w(i);
                    if(min_r > ratio) {
                        min_r = ratio;
                        min_i = m-n+i;
                    }
                }

                // 列の入れ替え
                for(int i=0; i<n; i++) {
                    std::swap(A(i, min_i), A(i, min_k));
                }
                std::swap(c(min_i), c(min_k));
                std::swap(idx[min_i], idx[min_k]);
            }

            xmin = Vector64f::zeros(m);
            for(int i=0; i<m; i++) {
                if(idx[i] < m-n) {
                    xmin(idx[i]) = x(i-m+n);
                }
            }
            fmin = c_.dot(xmin);

            delete[] idx;

        }
    }
}
