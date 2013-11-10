#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
using namespace std;

#include "../funopttoolkit/funopttoolkit.hpp"
using namespace funopt;

int main(int argc, char** argv) {
    srand((unsigned long)time(0));

    int d = 100;
    double elemA[] = {1.0, 2.0, 3.0,
                      2.0, 2.0, 3.0,
                      3.0, 3.0, 3.0}; 
    Matrix64f A(elemA, 3, 3);
    cout << A << endl;

    double elemb[] = {6.0, 7.0, 9.0};
    Vector64f b(elemb, 3);
    Vector64f x = A.solve(b, FUNOPT_FACTOR_LU);
    cout << x << endl;
}
