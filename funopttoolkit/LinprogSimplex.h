#ifndef _LINPROG_SIMPLEX_H_
#define _LINPROG_SIMPLEX_H_

#include "LinprogSolverBase.h"

namespace funopt {
    namespace linprog {
        class Simplex : public SolverBase {
        public:
            Simplex();
            virtual void solve(const Vector64f& c, const Matrix64f& A, const Vector64f& b);
        };
    }
}

#endif
