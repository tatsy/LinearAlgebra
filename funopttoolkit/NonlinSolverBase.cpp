#include "NonlinSolverBase.h"
#include "Vector64f.h"

namespace funopt {
    namespace nonlin {
        SolverBase::SolverBase() :
            f_ptr(0)
        {
        }

        SolverBase::SolverBase(funcNd* func_ptr) :
            f_ptr(func_ptr)
        {
        }

        SolverBase::SolverBase(const SolverBase& sb) :
            f_ptr(sb.f_ptr)
        {
        }

        SolverBase::~SolverBase()
        {
        }

        SolverBase& SolverBase::operator=(const SolverBase& sb)
        {
            f_ptr = sb.f_ptr;
            return *this;
        }

        Vector64f SolverBase::solve(const Vector64f& x0, const int maxiter, const double tol)
        {
            Vector64f x_opt;
            solve(x0, x_opt, maxiter, tol);
            return x_opt;
        }
    }
}
