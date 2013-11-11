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

	cout << "Input" << endl;
    cout << A << endl << endl;

	cout << "det(A) = " << A.det() << endl << endl;

	Matrix64f invA = A.inv();
	cout << "Inverse" << endl;
	cout << invA << endl << endl;

	cout << (A * invA) << endl;

    double elemb[] = {6.0, 7.0, 9.0};
    Vector64f b(elemb, 3);
    Matrix64f x = A.solve((Matrix64f)b, FUNOPT_FACTOR_QR);
    cout << x << endl;
}
