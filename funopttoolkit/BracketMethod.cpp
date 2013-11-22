#include <iostream>
using namespace std;

#include "BracketMethod.h"

static double sgn(double d)
{
    if(d < -1.0e-20) return -1.0;
    if(d >  1.0e-20) return  1.0;
    return 0.0;
}

namespace funopt {
    namespace nonlin {
        const double BracketMethod::GOLD   = 1.618034;
        const double BracketMethod::GLIMIT = 100.0;
        const double BracketMethod::TINY   = 1.0e-20;

        BracketMethod::BracketMethod() :
            lower(0.0),
            mid(0.0),
            upper(0.0)
        {
        }

        inline double BracketMethod::get_lower() const
        {
            return lower;
        }

        inline double BracketMethod::get_mid() const
        {
            return mid;
        }

        inline double BracketMethod::get_upper() const
        {
            return upper;
        }

        void BracketMethod::bracket(const double a, const double b, func1d* f_ptr)
        {
            double ax = a;
            double bx = b;
            double fa = (*f_ptr)(ax);
            double fb = (*f_ptr)(bx);
            if(fb > fa) {
                std::swap(ax, bx);
                std::swap(fa, fb);
            }

            double cx = b + GOLD * (bx - ax);
            double fc = (*f_ptr)(cx);
            while(fb > fc) {
                // ‘o‹Èü•âŠÔ
                double r = (bx - ax) * (fb - fc);
                double q = (bx - cx) * (fb - fa);
                double u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * sgn(q -r) * max(abs(q-r), TINY));
                double ulim = bx + GLIMIT * (cx - bx);
                double fu = 0.0;

                // ’l‚ÌXV
                if((bx - u) * (u - cx) > 0.0) {
                    fu = (*f_ptr)(u);
                    if(fu < fc) {
                        ax = bx;
                        bx = u;
                        fa = fb;
                        fb = fu;
                        break;
                    }
                    else if(fu > fb) {
                        cx = u;
                        fc = fu;
                        break;
                    }
                    u  = cx + GOLD * (cx - bx);
                    fu = (*f_ptr)(u);
                }
                else if((cx - u) * (u - ulim) > 0.0) {
                    fu = (*f_ptr)(u);
                    if(fu < fc) {
                        shift3(bx, cx, u, u + GOLD * (u - cx));
                        shift3(fb, fc, fu, (*f_ptr)(u));
                    }
                }
                else if((u - ulim) * (ulim - cx) >= 0.0) {
                    u  = ulim;
                    fu = (*f_ptr)(u);
                }
                else {
                    u  = cx + GOLD * (cx - bx);
                    fu = (*f_ptr)(u);
                }
                shift3(ax, bx, cx, u);
                shift3(fa, fb, fc, fu); 
            }

            lower = ax;
            mid   = bx;
            upper = cx;
        }

        inline void BracketMethod::shift2(double &a, double &b, const double c)
        {
            a = b;
            b = c;
        }

        inline void BracketMethod::shift3(double &a, double &b, double &c, const double d)
        {
            a = b;
            b = c;
            c = d;
        }

        inline void BracketMethod::mov3(double &a, double &b, double &c, const double d, const double e, const double f)
        {
            a = d;
            b = e;
            c = f;
        }
    }
}
