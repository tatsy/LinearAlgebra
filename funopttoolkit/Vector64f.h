#ifndef _VECTOR_64F_
#define _VECTOR_64F_

#include <iostream>
#include <sstream>
#include <cstring>
using namespace std;

#include "dll_macros.h"

namespace funopt {
    class Matrix64f;

    class DLL_EXPORT Vector64f {
        friend class Matrix64f;
        friend ostream& operator<<(ostream& os, const Vector64f& v);

    private:
        int ndim;
        double* data;

    public:
        Vector64f();
        Vector64f(int dim);
        Vector64f(const Vector64f& v);
        ~Vector64f();

        Vector64f& operator=(const Vector64f& v);

        double& operator()(int i);
        double operator()(int i) const;
        
        Vector64f operator+(const Vector64f& v) const;
        Vector64f operator-(const Vector64f& v) const;

        int dim() const;
        double norm() const;
        double norm2() const;
    };

    inline ostream& operator<<(ostream& os, const Vector64f& v)
    {
        os << "[ ";
        int lim = 10 < v.dim() ? 10 : v.dim();
        for(int i=0; i<lim; i++) {
            os << v(i);
            if(i != lim-1) os << ", ";
        }

        if(lim != v.dim()) {
            os << " ... and more";
        }
        os << " ]";
        return os;
    }

}


#endif
