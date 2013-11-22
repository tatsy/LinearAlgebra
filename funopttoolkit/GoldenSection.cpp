#include <cmath>

#include "GoldenSection.h"

namespace funopt {
    namespace nonlin {
        GoldenSection::GoldenSection(const double tol_) :
            xmin(0.0),
            fmin(0.0),
            tol(tol_)
        {
        }

        double GoldenSection::minimize(const func1d& func)
        {
            static const double R = 0.61803399;
            static const double C = 1.0 - R;

            double ax = get_lower();
            double bx = get_mid();
            double cx = get_upper();

            double x1, x2;
            double x0 = ax;
            double x3 = cx;
            if(abs(cx - bx) > abs(bx - ax)) {
                x1 = bx;
                x2 = bx + C * (cx - bx);
            } else {
                x2 = bx;
                x1 = bx - C * (bx - ax);
            }
            double f1 = func(x1);
            double f2 = func(x2);
            while(abs(x3 - x0) > tol * (abs(x1) + abs(x2))) {
                if(f2 < f1) {
                    shift3(x0, x1, x2, R * x2 + C * x3);
                    shift2(f1, f2, func(x2));
                } else {
                    shift3(x3, x2, x1, R * x1 + C * x0);
                    shift2(f2, f1, func(x1));
                }
            }

            if(f1 < f2) {
                xmin = x1;
                fmin = f1;
            } else {
                xmin = x2;
                fmin = f2;
            }
            return xmin;
        }
    }
}

