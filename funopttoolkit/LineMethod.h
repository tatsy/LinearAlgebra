#ifndef _LINE_METHOD_H_
#define _LINE_METHOD_H_

#include "Vector64f.h"
#include "NonlinFunc.h"
#include "NonlinSolverBase.h"

namespace funopt {
    namespace nonlin {
        class LineMethod : public SolverBase {
        protected:
            double fmin;
            Vector64f xmin;

        public:
            LineMethod();
            double linmin(const funcNd& func, const Vector64f& x, const Vector64f& e, const int maxiter, const double tol);

            virtual Vector64f& get_xmin();
            virtual double     get_fmin() const;
        };
    }
}

#endif
