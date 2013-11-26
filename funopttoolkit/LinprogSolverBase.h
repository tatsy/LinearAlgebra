#ifndef _LINPROG_SOLVER_BASE_H_
#define _LINPROG_SOLVER_BASE_H_

#include "Vector64f.h"

namespace funopt {
    class Matrix64f;

    namespace linprog {
        class SolverBase {
        protected:
            double    fmin;
            Vector64f xmin;

        public:
            SolverBase();

            virtual const Vector64f& getXMin() const;
            virtual double getFMin() const;

            // solve linear programming such as,
            // minimize: (c,x)
            // subject to: Ax = b
            virtual void solve(const Vector64f& c, const Matrix64f& A, const Vector64f& b) = 0; 
        };
    }
}

#endif
