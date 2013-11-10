#ifndef _LINSOLVE_H_
#define _LINSOLVE_H_

#include "dll_macros.h"

namespace funopt {
    class Vector64f;
    class Matrix64f;

    void factor_lu(Matrix64f& A, Matrix64f& LU, int* order);
    void solve_lu(Matrix64f& A, Vector64f& b, Vector64f& x);
}

#endif
