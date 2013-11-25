#ifndef _NONLIN_CONJUGATE_GRADIENT_H_
#define _NONLIN_CONJUGATE_GRADIENT_H_

#include "LineMethod.h"

namespace funopt {
    namespace nonlin {
        class ConjugateGradient : public LineMethod {
        public:
            ConjugateGradient();
            void solve(const funcNd& func, const Vector64f& x0, Vector64f& x_opt, const int maxiter, const double tol);
        };
    }
}

#endif
