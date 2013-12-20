#include <iostream>

#include "ConjugateGradient.h"

namespace funopt {
    namespace nonlin {
        ConjugateGradient::ConjugateGradient() 
        {
        }

        void ConjugateGradient::solve(const funcNd& func, const Vector64f& x0, Vector64f& x_opt, const int maxiter, const double tol)
        {
            static const double EPS  = 1.0e-18;
            static const double GTOL = 1.0e-8;

            const int n = x0.dim();
            Vector64f x = x0;

            double fx = func(x);
            Vector64f e = -func.grad(x);
            Vector64f g = e;
            Vector64f h = e;

            for(int it=0; it<maxiter; it++) {
                double fret = linmin(func, x, e, maxiter, tol);
                if(2.0 * abs(fret - fx) <= tol * (abs(fret) * abs(fx) + EPS)) {
                    break;
                }

                fx = fret;
                x  = xmin;

                g = func.grad(x);
                double test = 0.0;
                double den = std::max(abs(fx), 1.0);
                for(int j=0; j<n; j++) {
                    double temp = abs(g(j)) * std::max(abs(x(j)), 1.0) / den;
                    if(temp > test) test = temp;
                }

                if(test < GTOL) {
                    break;
                }

                double gg = e.dot(e);
                double dgg = (g + e).dot(e);
                if(gg == 0.0) {
                    break;
                }

                double gam = dgg / gg;
                e = -g;
                g = h = e + gam * h;
            }
            x_opt = x;
        }
    }
}
