#include <cstdio>
#include <cstdlib>
#include <ctime>

#include "../funopttoolkit/funopttoolkit.hpp"
using namespace funopt;

int main(int argc, char** argv) {
    srand((unsigned long)time(0));

    int d = 1000;
    Matrix64f A(d, d);
    for(int i=0; i<d; i++) {
        for(int j=0; j<d; j++) {
            A(i, j) = rand() % 1000;
        }
    }

    Vector64f x(d);
    for(int i=0; i<d; i++) {
        x(i) = rand() % 1000;
    }
    Vector64f b = A * x;
   

    Vector64f y = A.solve(b, FUNOPT_FACTOR_LU);
    printf("%f\n", (x - y).norm());
}
