#include "Vector64f.h"
#include "Matrix64f.h"
#include "LinprogSolverBase.h"

namespace funopt {
    namespace linprog {
        SolverBase::SolverBase() :
            fmin(0.0),
            xmin()
        {
        }

        const Vector64f& SolverBase::getXMin() const
        {
            return xmin;
        }

        double SolverBase::getFMin() const
        {
            return fmin;
        }
    }
}
    