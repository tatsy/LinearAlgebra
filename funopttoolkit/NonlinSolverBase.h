#ifndef _NONLIN_SOLVER_BASE_H_
#define _NONLIN_SOLVER_BASE_H_

#include "funopt_macros.h"

namespace funopt {
    class Vector64f;

    namespace nonlin {
        // ソルバの種類
        enum SolverType {
            SOLVER_NELDER_SIMPLEX
        };
        
        // 関数クラス (継承して使う)
        class func1d {
        public:
            // 関数の評価
            virtual double operator()(double x) = 0;

            // 関数の微分
            virtual double df(double x) {
                massert(false, "derivative not defined");
            }
        };

        // 関数クラス (継承して使う)
        class funcNd {
        public:
            // 関数の評価
            virtual double operator()(const Vector64f& x) = 0;
            // 関数の勾配
            virtual double grad(const Vector64f& x) {
                massert(false, "gradient not defined");
            }
        };

        // ソルバ基底クラス
        class SolverBase {
        protected:
            funcNd* f_ptr;

        public:
            SolverBase();
            explicit SolverBase(funcNd* func_ptr);
            SolverBase(const SolverBase& sb);
            ~SolverBase();

            SolverBase& operator=(const SolverBase& sb);

            Vector64f solve(const Vector64f& x0, const int maxiter, const double tol);
            virtual void solve(const Vector64f& x0, Vector64f& x_opt, const int maxiter, const double tol) = 0;
        };
    };
};

#endif
