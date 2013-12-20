#ifndef _MATRIX_64F_
#define _MATRIX_64F_

#ifdef __MAT64F_EXPORT__
#define MAT64F_DLL_EXPORT __declspec(dllexport)
#else
#define MAT64F_DLL_EXPORT __declspec(dllimport)
#endif

#include <iostream>
#include <vector>
#include <sstream>
#include <cstring>

#include "funopt_enum.h"
#include "matrix_enums.h"

namespace funopt {
    class Vector64f;

    class MAT64F_DLL_EXPORT Matrix64f {
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

		// 乱数行列
		static Matrix64f rand(int rows, int cols);

        // 対角行列
        static Matrix64f diag(std::vector<double>& v);

        // デストラクタ
        ~Matrix64f();

        Matrix64f& operator=(const Matrix64f& m);

        double& operator()(int i, int j);
        double  operator()(int i, int j) const;
        Matrix64f& operator+=(const Matrix64f& A);
        Matrix64f& operator-=(const Matrix64f& A);
        
        int rows() const;
        int cols() const;

        Vector64f row(int i) const;
        Vector64f col(int j) const;

        // 部分行列の取り出し
        Matrix64f submat(int i, int j, int rows, int cols) const;

		// 行列の転置
		Matrix64f trans() const;

		// 行列式を求める
		double det() const;

		// 行列ノルム
        double norm(MatrixNormType type = MAT_NORM_FROBENIUS) const;

		// 逆行列を求める
		Matrix64f inv() const;

        // 線形問題を解く
		Matrix64f solve(const Matrix64f& b, SolverType solver_type=SOLVER_LU);

        // 固有値を求める
        void eig(Matrix64f& val, Matrix64f& vec) const;

	private:
		// LU分解
		void factor_lu(Matrix64f& LU, int* order) const;

		// LU分解により線形問題を解く
		void solve_lu(const Matrix64f& b, Matrix64f& x) const;

		// QR分解
		void factor_qr(Matrix64f& Q, Matrix64f& R) const;

		// QR分解により線形問題を解く
		void solve_qr(const Matrix64f& b, Matrix64f& x) const;

        // CG法により線形問題を解く
        void solve_cg(const Matrix64f& b, Matrix64f& x) const;
    };
}

MAT64F_DLL_EXPORT std::ostream& operator<<(std::ostream& os, const funopt::Matrix64f& m);


// 演算子定義

MAT64F_DLL_EXPORT funopt::Matrix64f operator+(const funopt::Matrix64f& A, const funopt::Matrix64f& B);
MAT64F_DLL_EXPORT funopt::Matrix64f operator-(const funopt::Matrix64f& A, const funopt::Matrix64f& B);
MAT64F_DLL_EXPORT funopt::Matrix64f operator*(const funopt::Matrix64f& A, const funopt::Matrix64f& B);
MAT64F_DLL_EXPORT funopt::Vector64f operator*(const funopt::Matrix64f& A, const funopt::Vector64f& v);
MAT64F_DLL_EXPORT funopt::Matrix64f operator*(const funopt::Matrix64f& A, double d);
MAT64F_DLL_EXPORT funopt::Matrix64f operator*(double d, const funopt::Matrix64f& A);
MAT64F_DLL_EXPORT funopt::Matrix64f operator/(const funopt::Matrix64f& A, double d);

#endif
