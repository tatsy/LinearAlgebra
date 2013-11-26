#ifndef _LINPROG_SOLVER_H_
#define _LINPROG_SOLVER_H_

#ifndef __LINPROG_EXPORT__
#define LINPROG_EXPORT __declspec(dllexport)
#else
#define LINPROG_EXPORT __declspec(dllimport)
#endif

#include "Vector64f.h"
#include "LinprogSolverBase.h"

namespace funopt {
    class Matrix64f;

    namespace linprog {
        enum SolverType {
            SOLVER_SIMPLEX
        };

        class LINPROG_EXPORT Solver {
        private:
            SolverBase* base;

        public:
            Solver();
            ~Solver();

            const Vector64f& getXMin() const;
            double getFMin() const;

            void solve(const Vector64f& c, const Matrix64f& A, const Vector64f& b, SolverType type);

        private:
            Solver(const Solver& solver);
            Solver& operator=(const Solver& solver);
        };
    }
}

#endif

