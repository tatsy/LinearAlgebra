#include "Vector64f.h"
#include "Matrix64f.h"
#include "QuasiNewton.h"

namespace funopt {
    namespace nonlin {
        QuasiNewton::QuasiNewton()
        {
        }

        void QuasiNewton::solve(const funcNd& func, const Vector64f& x0, Vector64f& x_opt, const int maxiter, const double tol)
        {
            static const double EPS   = 1.0e-20;
            static const double TOLX  = EPS * 4.0;
            static const double STPMX = 0.01;

            double fret;
            bool check;
            const int n = x0.dim();

            Vector64f x, xnew;
            x = xnew = x0;

            Matrix64f hessian = Matrix64f::eye(n);
            double fx = func(x);
            Vector64f g = func.grad(x);
            Vector64f e = -g;
            double sum = x.norm();
            double stpmax = STPMX * max(sum, (double)n);
            for(int it=0; it<maxiter; it++) {
                linsearch(x, fx, g, e, xnew, fret, stpmax, check, func);
                fx = fret;
                e  = xnew - x;
                x  = xnew;

                double test = 0.0;
                for(int i=0; i<n; i++) {
                    double temp = abs(e(i)) / max(abs(x(i)), 1.0);
                    if(temp > test) {
                        test = temp;
                    }
                }

                if(test < TOLX) {
                    break;
                }

                Vector64f dg = g;
                g = func.grad(x);

                test = 0.0;
                double den = max(abs(fret), 1.0);
                for(int i=0; i<n; i++) {
                    double temp = abs(g(i)) * max(abs(x(i)), 1.0) / den;
                    if(temp > test) test = temp;
                }

                if(test < tol) {
                    break;
                }

                dg = g - dg;
                Vector64f hdg = hessian * dg;

                double fad = 0.0;
                double fac = dg.dot(e);
                double fae = dg.dot(hdg);
                double sumdg = dg.dot(dg);
                double sume = e.dot(e);
                if(fac > sqrt(EPS * sumdg * sume)) {
                    fac = 1.0 / fac;
                    fad = 1.0 / fae;
                    dg  = fac * e - fad * hdg;
                    hessian += fac * e.tensor(e) - fad * hdg.tensor(hdg) + fae * dg.tensor(dg);
                }
                e = -(hessian * g);
            }
            x_opt = x;
        }

        void QuasiNewton::linsearch(Vector64f& xold, const double fold, Vector64f& g, Vector64f& p, Vector64f& x, double& f, const double stpmax, bool & check, const funcNd& func)
        {
            static const double ALF  = 1.0e-4;
            static const double TOLX = 1.0e-20;
            const int n= xold.dim();
            check = false;
            double sum = p.norm();

            if(sum > stpmax) {
                p *= stpmax / sum;
            }
            double slope = g.dot(p);
            //massert(slope < 0.0, "Roundoff problem in linsearch");

            double test = 0.0;
            for(int i=0; i<n; i++) {
                double temp = abs(p(i)) / max(abs(xold(i)), 1.0);
                if(temp > test) {
                    test = temp;
                }
            }

            double tmplam = 0.0;
            double alam2 = 0.0, f2 = 0.0;
            double alamin = TOLX / test;
            double alam   = 1.0;
            for(;;) {
                x = xold + alam * p;
                f = func(x);
                if(alam < alamin) {
                    x = xold;
                    check = true;
                    return;
                } else if(f <= fold + ALF * alam * slope) {
                    return;
                } else {
                    if(alam == 1.0) {
                        tmplam = -slope / (2.0 * (f - fold - slope));
                    } else {
                        double rhs1 = f - fold - alam * slope;
                        double rhs2 = f2 - fold - alam2 * slope;
                        double a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
                        double b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
                        if(a == 0.0) {
                            tmplam = -slope / (2.0 * b);
                        } else {
                            double disc = b * b - 3.0 * a * slope;
                            if(disc < 0.0) {
                                tmplam = 0.5 * alam;
                            } else if(b <= 0.0) {
                                tmplam = (-b + sqrt(disc)) / (3.0 * a);
                            } else {
                                tmplam = -slope / (b + sqrt(disc));
                            }
                        }
                        if(tmplam > 0.5 * alam) {
                            tmplam = 0.5 * alam;
                        }
                    }
                }
                alam2 = alam;
                f2    = f;
                alam  = max(tmplam, 0.1 * alam);
            }
        }
    }
}
