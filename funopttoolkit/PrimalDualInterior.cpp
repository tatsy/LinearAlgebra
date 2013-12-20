#include <iostream>

#include "PrimalDualInterior.h"
#include "Matrix64f.h"

namespace funopt {
    namespace linprog {
        PrimalDualInterior::PrimalDualInterior() :
            SolverBase()
        {
        }

        void PrimalDualInterior::solve(const Vector64f& c_, const Matrix64f& A_, const Vector64f& b_)
        {
            Vector64f c  = c_;
            Matrix64f A  = A_;
            Matrix64f At = A.trans();
            Vector64f b  = b_;

            static const int    maxiter = 200;
            static const double EPS     = 1.0e-6;
            static const double SIGMA   = 0.9;
            static const double DELTA   = 0.02;
            static const double BIG     = 1.0e20;

            double primal_obj, dual_obj;

            int status = -1;
            int m = A.rows();
            int n = A.cols();

            double rpfact = 1.0 + b.norm();
            double rdfact = 1.0 + c.norm();
            Vector64f oneN = Vector64f::ones(n);
            Vector64f x = 1000.0 * oneN;
            Vector64f z = 1000.0 * oneN;
            Vector64f y = 1000.0 * Vector64f::ones(m);

            double normrp_old = BIG;
            double normrd_old = BIG;

            for(int it=0; it<maxiter; it++) {
                Vector64f ax = A * x;
                Vector64f rp = ax - b;
                double normrp = rp.norm() / rpfact;
                Vector64f aty = At * y;
                Vector64f rd = aty + z - c;
                double normrd = rd.norm() / rdfact;
                double gamma = x.dot(z);
                double mu = DELTA * gamma / n;
                primal_obj = c.dot(x);
                dual_obj   = b.dot(y);
                double gamma_norm = gamma / (1.0 + abs(primal_obj));
                if(normrp < EPS && normrd < EPS && gamma_norm < EPS) {
                    status = 0;
                    break;
                }

                if(normrp > 1000.0 * normrp_old && normrp > EPS) {
                    status = 1;
                    break;
                }

                if(normrd > 1000.0 * normrd_old && normrd > EPS) {
                    status = 2;
                    break;
                }

                Vector64f d = x / z;
                Matrix64f adat = A * d.asDiag() * At;

                Vector64f tempn = x - mu / z - d * rd;
                Vector64f tempm = A * tempn;
                Vector64f rhs = -rp + tempm;
                Vector64f dy = adat.inv() * rhs;
                tempn = At * dy;
                Vector64f dz = -tempn - rd;
                Vector64f dx = -d * dz + mu / z - x; 

                double alpha_p = 1.0;
                for(int j=0; j<n; j++) {
                    if(x(j) + alpha_p * dx(j) < 0.0) {
                        alpha_p = -x(j) / dx(j);
                    }
                }

                double alpha_d = 1.0;
                for(int j=0; j<n; j++) {
                    if(z(j) + alpha_d * dz(j) < 0.0) {
                        alpha_d = -z(j) / dz(j);
                    }
                }

                alpha_p = std::min(alpha_p * SIGMA, 1.0);
                alpha_d = std::min(alpha_d * SIGMA, 1.0);

                x += alpha_p * dx;
                z += alpha_d * dz;
                y += alpha_d * dy;

                normrp_old = normrp;
                normrd_old = normrd;
            }
            status = 3;
            fmin = primal_obj;
            xmin = x;
        }
    }
}
