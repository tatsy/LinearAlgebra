#include "NonlinSolverBase.h"
#include "Vector64f.h"

namespace funopt {
    namespace nonlin {
        SolverBase::SolverBase()
        {
        }

        Vector64f SolverBase::solve(const funcNd& func, const Vector64f& x0, const int maxiter, const double tol)
        {
            Vector64f x_opt;
            solve(func, x0, x_opt, maxiter, tol);
            return x_opt;
        }
    }
}
