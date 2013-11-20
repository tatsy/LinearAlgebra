#include <iostream>
using namespace std;

#include "../../funopttoolkit/funopttoolkit.hpp"
using namespace funopt;

int main(int argc, char** argv) {
	//double elemA[] = {1, 2, 3, 4,
	//                  2, 2, 3, 4,
 //                     3, 3, 3, 4,
 //                     4, 4, 4, 4 };
 //   Matrix64f A(elemA, 4, 4);

	 Matrix64f A = Matrix64f::rand(100, 100);
	 A = A + A.trans();
	
    Matrix64f L, U;
	A.eig(L, U);

	// cout << L << endl << endl;
	// cout << U << endl << endl;
}
