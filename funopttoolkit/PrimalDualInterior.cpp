#include "PrimalDualInterior.h"

#include "Matrix64f.h"

namespace funopt {
    namespace linprog {
        PrimalDualInterior::PrimalDualInterior() :
            SolverBase()
        {
        }

        void PrimalDualInterior::solve(const Vector64f& c_, const Matrix64f& A_, const Vector64f& b_)
        {
            Vector64f c = c_;
            Matrix64f A = A_;
            Vector64f b = b_;
        }
    }
}
