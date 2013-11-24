#ifndef _POWELL_H_
#define _POWELL_H_

#include "LineMethod.h"

namespace funopt {
    namespace nonlin {
        class Powell : public LineMethod {
        public:
            Powell();
            double minimize(const funcNd& func, const Vector64f& x0, const int maxiter, const double tol);
            void solve(const funcNd& func, const Vector64f& x0, Vector64f& x_opt, const int maxiter, const double tol);
        };
    }
}

#endif
