#include "../../funopttoolkit/funopttoolkit.hpp"
using namespace funopt;

int main(int argc, char** argv) {
    double elemA[] = { 1.0, 2.0, 3.0,
                       2.0, 2.0, 3.0,
                       3.0, 3.0, 3.0 };
    
    Matrix64f A(elemA, 3, 3);
    double fro = A.norm(MAT_NORM_FROBENIUS);
    double spe = A.norm(MAT_NORM_SPECTRAL);

    printf("frobenius: %f\n", fro);
    printf(" spectral: %f\n", spe);
}
