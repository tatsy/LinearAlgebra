#ifndef _NONLIN_SOLVER_BASE_H_
#define _NONLIN_SOLVER_BASE_H_

#include "funopt_macros.h"
#include "NonlinFunc.h"

namespace funopt {
    class Vector64f;

    namespace nonlin {
        // ソルバの種類
        enum SolverType {
            SOLVER_NELDER_SIMPLEX,
            SOLVER_POWELL,
            SOLVER_CONJUGATE_GRADIENT,
            SOLVER_QUASI_NEWTON
        };
        
        // ソルバ基底クラス
        class SolverBase {
        public:
            SolverBase();

            Vector64f solve(const funcNd& func, const Vector64f& x0, const int maxiter, const double tol);
            virtual void solve(const funcNd& func, const Vector64f& x0, Vector64f& x_opt, const int maxiter, const double tol) = 0;
        };
    };
};

#endif
