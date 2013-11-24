#include "LineMethod.h"
#include "Brent.h"

namespace funopt {
    namespace nonlin {
        LineMethod::LineMethod() :
            fmin(0.0),
            xmin()
        {
        }

        double LineMethod::linmin(const funcNd& func, const Vector64f& x, const Vector64f& e, const int maxiter, const double tol) 
        {
            const int n = x.dim();
            
            func1d f1dlin(func, x, e);
            double ax = 0.0;
            double xx = 1.0;

            Brent brent;
            brent.bracket(ax, xx, f1dlin);
            double xm = brent.minimize(f1dlin, maxiter, tol); 
            xmin = Vector64f(n);
            for(int j=0; j<n; j++) {
                xmin(j) = x(j) + xm * e(j);
            }
            fmin = brent.get_fmin();
            return fmin;
        }

        Vector64f& LineMethod::get_xmin()
        {
            return xmin;
        }

        double LineMethod::get_fmin() const
        {
            return fmin;
        }
    }
}
