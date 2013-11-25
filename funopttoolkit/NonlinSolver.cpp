#define __NONLIN_EXPORT__
#include "NonlinSolver.h"
#include "NelderSimplex.h"
#include "Powell.h"
#include "ConjugateGradient.h"
#include "Vector64f.h"

namespace funopt {
    namespace nonlin {
        Solver::Solver() 
        {
        }

        void Solver::solve(const funcNd& func, const Vector64f& x0, Vector64f& x_opt, SolverType type, const int maxiter, const double tol)
        {
            SolverBase* base = 0;
            switch(type) {
            case SOLVER_NELDER_SIMPLEX:
                base = new NelderSimplex();
                break;

            case SOLVER_POWELL:
                base = new Powell();
                break;
            
            case SOLVER_CONJUGATE_GRADIENT:
                base = new ConjugateGradient();
                break;
            }
            base->solve(func, x0, x_opt, maxiter, tol);
        }

        Vector64f Solver::solve(const funcNd& func, const Vector64f& x0, SolverType type, const int maxiter, const double tol)
        {
            Vector64f x_opt;
            solve(func, x0, x_opt, type, maxiter, tol);
            return x_opt;
        }
    }
}
