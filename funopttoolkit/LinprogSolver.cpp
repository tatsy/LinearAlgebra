#include "Vector64f.h"
#include "Matrix64f.h"

#include "LinprogSimplex.h"
#include "PrimalDualInterior.h"

#define __LINPROG_EXPORT__
#include "LinprogSolver.h"

namespace funopt {
    namespace linprog {
        Solver::Solver() :
            base(0)
        {
        }

        Solver::~Solver()
        {
            delete base;
        }

        const Vector64f& Solver::getXMin() const
        {
            return base->getXMin();
        }

        double Solver::getFMin() const
        {
            return base->getFMin();
        }

        void Solver::solve(const Vector64f& c, const Matrix64f& A, const Vector64f& b, SolverType type)
        {
            switch(type) {
            case SOLVER_SIMPLEX:
                base = new Simplex();
                break;

            case SOLVER_PRIMAL_DUAL:
                base = new PrimalDualInterior();
                break;
            }
            base->solve(c, A, b);
        }
    }
}
