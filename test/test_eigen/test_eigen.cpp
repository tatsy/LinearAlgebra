#include <iostream>
using namespace std;

#include "../../funopttoolkit/funopttoolkit.hpp"
using namespace funopt;

int main(int argc, char** argv) {
    /*
    double elemA[] = {
        1, 2, 3, 4,
        2, 2, 3, 4,
        3, 3, 3, 4,
        4, 4, 4, 4 
    };

    Matrix64f A(elemA, 4, 4);
    */

    const int n = 15;
    Matrix64f A = Matrix64f::rand(n, n);
    A = A + A.trans();
    Matrix64f L, U;
    A.eig(L, U);

    Matrix64f B = U * L * U.trans();
    double err = 0.0;
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            double d = A(i, j) - B(i, j);
            err += d * d;
        }
    }
    printf("error = %.15f\n", err);
    
    cout << L << endl << endl;
    //cout << U << endl << endl;
    //cout << B << endl << endl;
    //cout << A << endl << endl;
}
