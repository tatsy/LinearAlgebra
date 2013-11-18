#include <iostream>
using namespace std;

#include "../../funopttoolkit/funopttoolkit.hpp"
using namespace funopt;

int main(int argc, char** argv) {
	Matrix64f A = Matrix64f::rand(5, 5);
	A = A + A.trans();
	Matrix64f L, U;
	A.eig(L, U);

	cout << L << endl << endl;
	cout << U << endl << endl;
}
