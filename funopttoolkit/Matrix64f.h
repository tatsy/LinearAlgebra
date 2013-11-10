#ifndef _MATRIX_64F_
#define _MATRIX_64F_

#include <iostream>
#include <sstream>
#include <cstring>
using namespace std;

#include "dll_macros.h"

namespace funopt {
    class Vector64f;

    class DLL_EXPORT Matrix64f {
        friend ostream& operator<<(ostream& os, const Matrix64f& m);


    private:
        int nrows, ncols;
        double* data;

    public:
        // コンストラクタ
        Matrix64f();
        Matrix64f(int rows, int cols);
        Matrix64f(double* data, int rows, int cols);
        Matrix64f(const Matrix64f& m);

        // デストラクタ
        ~Matrix64f();

        Matrix64f& operator=(const Matrix64f& m);

        double& operator()(int i, int j);
        double  operator()(int i, int j) const;

        Vector64f operator*(const Vector64f& v);

        int rows() const;
        int cols() const;

        // 線形問題を解く
        Vector64f solve(Vector64f& b, int decomp_type);    
    };

    inline ostream& operator<<(ostream& os, const Matrix64f& m)
    {
        for(int i=0; i<m.rows(); i++) {
            os << (i == 0 ? "[" : " ");
            os << "[ ";
            for(int j=0; j<m.cols(); j++) {
                os << m(i, j) << (j == m.cols()-1 ? " ]" : ", ");
            }
            os << (i == m.rows()-1 ? "]" : "\n");
        }
        return os;
    }
}

#endif
