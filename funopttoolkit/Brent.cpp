#include <cmath>

#define __BRENT_EXPORT__
#include "Brent.h"

namespace funopt {
    namespace nonlin {
        Brent::Brent() :
            xmin(0.0),
            fmin(0.0)
        {
        }

        double Brent::minimize(const func1d& func, const int maxiter, const double tol)
        {
            static const double CGOLD = 0.3819660;
            static const double ZEPS  = 1.0e-12;
            double d, e;
            d = e = 0.0;

            double ax = get_lower();
            double bx = get_mid();
            double cx = get_upper();

            double a, b;
            a = ax < cx ? ax : cx;
            b = ax > cx ? ax : cx;

            double x, w, v;
            x = w = v = bx;

            double fx, fw, fv;
            fx = fw = fv = func(x);

            for(int it=0; it<maxiter; it++) {
                double xm = 0.5 * (a + b);
                double tol1 = tol * abs(x) + ZEPS;
                double tol2 = 2.0 * tol1;
                
                // I—¹”»’è
                if(abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                    break;
                }

                // ‘o‹Èü•âŠÔ
                if(abs(e) > tol1) {
                    double r = (x - w) * (fx - fv);
                    double q = (x - v) * (fx - fw);
                    double p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);
                    if(q > 0.0) {
                        p = -p;
                    }
                    q = abs(q);
                    double etemp = e;
                    e = d;
                    if(abs(p) >= abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                        e = x >= xm ? a - x : b - x;
                        d = CGOLD * e;
                    }
                    else {
                        d = p / q;
                        double u = x + d;
                        if(u - a < tol2 || b - u < tol2) {
                            d = sgn(xm - x) * tol1;
                        }
                    }
                }
                else {
                    e = x >= xm ? a - x : b - x;
                    d = CGOLD * e;
                }
                double u = abs(d) >= tol1 ? x + d : x + sgn(d) * tol1;
                double fu = func(u);

                if(fu <= fx) {
                    if(u >= x) {
                        a = x;
                    } else {
                        b = x;
                    }
                    shift3(v, w, x, u);
                    shift3(fv, fw, fx, fu);
                }
                else {
                    if(u < x) {
                        a = u;
                    } else {
                        b = u;
                    }
                    if(fu <= fw || w == x) {
                        v = w;
                        w = u;
                        fv = fw; 
                        fw = fu;
                    } else if(fu <= fv || v == x || v == w) {
                        v = u;
                        fv = fu;
                    }
                }
            }

            fmin = fx;
            xmin = x;
            return xmin;        
        }

        double Brent::get_fmin() const
        {
            return fmin;
        }

        double Brent::get_xmin() const
        {
            return xmin;
        }
    }
}
