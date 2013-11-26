#ifndef _PRIMAL_DUAL_INTERIOR_H_
#define _PRIMAL_DUAL_INTERIOR_H_

#include "LinprogSolverBase.h"

namespace funopt {
    namespace linprog {
        class PrimalDualInterior : public SolverBase {
        public:
            PrimalDualInterior();
            virtual void solve(const Vector64f& c, const Matrix64f& A, const Vector64f& b);
        };
    }
}

#endif
