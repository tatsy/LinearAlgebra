#ifndef _MATRIX_64F_
#define _MATRIX_64F_

#include <iostream>
#include <sstream>
#include <cstring>
using namespace std;

#include "funopt_enum.h"
#include "dll_macros.h"
#include "matrix_enums.h"

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
		Matrix64f(const Vector64f& v);

		// 単位行列
		static Matrix64f eye(int n);

		// 零行列
		static Matrix64f zeros(int rows, int cols);

        // デストラクタ
        ~Matrix64f();

        Matrix64f& operator=(const Matrix64f& m);

        double& operator()(int i, int j);
        double  operator()(int i, int j) const;

		Matrix64f operator+(const Matrix64f& m) const;
		Matrix64f operator-(const Matrix64f& m) const;

        Vector64f operator*(const Vector64f& v) const;
		Matrix64f operator*(const Matrix64f& m) const;
		Matrix64f operator*(const double d) const;
		Matrix64f operator/(const double d) const;

        int rows() const;
        int cols() const;

		// 行列の転置
		Matrix64f trans() const;

		// 行列式を求める
		double det() const;

		// 行列ノルム
        double norm(MatrixNormType type = MAT_NORM_FROBENIUS) const;

		// 逆行列を求める
		Matrix64f inv() const;

        // 線形問題を解く
		Matrix64f solve(Matrix64f& b, int decomp_type=FUNOPT_FACTOR_LU);

        // 固有値を求める
        void eig(Matrix64f& eval, Matrix64f& evec) const;

	private:
		// LU分解
		void factor_lu(Matrix64f& LU, int* order) const;

		// LU分解により線形問題を解く
		void solve_lu(Matrix64f& b, Matrix64f& x) const;

		// QR分解
		void factor_qr(Matrix64f& Q, Matrix64f& R) const;

		// QR分解により線形問題を解く
		void solve_qr(Matrix64f& b, Matrix64f& x) const;
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
