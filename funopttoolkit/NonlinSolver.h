#ifndef _NONLIN_SOLVER_H_
#define _NONLIN_SOLVER_H_

#ifdef __NONLIN_EXPORT__
#define NONLIN_EXPORT __declspec(dllexport)
#else
#define NONLIN_EXPORT __declspec(dllimport)
#endif

#include "NonlinSolverBase.h"

namespace funopt {
    namespace nonlin {
        // Downhill Simplex Method
        class NONLIN_EXPORT Solver {
        public:
            Solver();
            Solver(const Solver& s);
            ~Solver();

            Solver& operator=(const Solver& s);

            Vector64f solve(funcNd* func_ptr, const Vector64f& x0, SolverType type, const int maxiter, const double tol);
            void solve(funcNd* func_ptr, const Vector64f& x0, Vector64f& x_opt, SolverType type, const int maxiter, const double tol);
        };
    }
}

#endif
