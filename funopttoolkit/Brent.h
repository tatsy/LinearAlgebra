#ifndef _BRENT_H_
#define _BRENT_H_

#include "BracketMethod.h"

#ifdef __BRENT_EXPORT__
#define BRENT_EXPORT __declspec(dllexport)
#else
#define BRENT_EXPORT __declspec(dllimport)
#endif

namespace funopt {
    namespace nonlin {
        class BRENT_EXPORT Brent : public BracketMethod {
        private:
            double xmin;
            double fmin;

        public:
            Brent();
            double minimize(const func1d& func, const int maxiter, const double tol);
            double get_fmin() const;
            double get_xmin() const;
        };
    }
}

#endif
