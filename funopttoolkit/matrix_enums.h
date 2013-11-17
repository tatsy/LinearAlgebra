#ifndef _MATRIX_ENUMS_H_
#define _MATRIX_ENUMS_H_

enum MatrixNormType {
    MAT_NORM_FROBENIUS,
    MAT_NORM_SPECTRAL
};

enum LinearSolverType {
    DECOMP_LU,
    DECOMP_QR,
    DECOMP_SVD,
    ITERATE_CG,
    ITERATE_BICG
};

#endif
