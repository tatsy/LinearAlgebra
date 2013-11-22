#ifndef _GOLDEN_SECTION_H_
#define _GOLDEN_SECTION_H_

#include "BracketMethod.h"

namespace funopt {
    namespace nonlin {
        class GoldenSection : public BracketMethod {
        private:
            double xmin, fmin;
            const double tol;

        public:
            GoldenSection(const double tol=3.0-8);
            double minimize(const func1d& func);
        };
    }
}

#endif
