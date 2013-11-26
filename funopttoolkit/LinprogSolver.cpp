#include "Vector64f.h"
#include "Matrix64f.h"

#include "LinprogSimplex.h"

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
            }
            base->solve(c, A, b);
        }
    }
}
