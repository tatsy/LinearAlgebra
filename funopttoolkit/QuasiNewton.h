#ifndef _QUASI_NEWTON_H_
#define _QUASI_NEWTON_H_

#include "NonlinSolverBase.h"

namespace funopt {
    namespace nonlin {
        class QuasiNewton : public SolverBase {
        public:
            QuasiNewton();
            void solve(const funcNd& func, const Vector64f& x0, Vector64f& x_opt, const int maxiter, const double tol);

        private:
            void linsearch(Vector64f& xold, const double fold, Vector64f& g, Vector64f& p, Vector64f& x, double& f, const double stpmax, bool & check, const funcNd& func);
        };
    }
}

#endif
