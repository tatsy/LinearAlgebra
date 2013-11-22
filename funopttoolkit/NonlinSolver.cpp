#define __NONLIN_EXPORT__
#include "NonlinSolver.h"
#include "NelderSimplex.h"
#include "Vector64f.h"

namespace funopt {
    namespace nonlin {
        Solver::Solver() 
        {
        }

        Solver::Solver(const Solver& s)
        {
        }

        Solver::~Solver()
        {
        }

        Solver& Solver::operator=(const Solver& s)
        {
            return *this;
        }

        void Solver::solve(funcNd* func_ptr, const Vector64f& x0, Vector64f& x_opt, SolverType type, const int maxiter, const double tol)
        {
            SolverBase* base = 0;
            switch(type) {
            case SOLVER_NELDER_SIMPLEX:
                base = new NelderSimplex(func_ptr);
            }
            base->solve(x0, x_opt, maxiter, tol);
        }

        Vector64f Solver::solve(funcNd* func_ptr, const Vector64f& x0, SolverType type, const int maxiter, const double tol)
        {
            Vector64f x_opt;
            solve(func_ptr, x0, x_opt, type, maxiter, tol);
            return x_opt;
        }
    }
}
