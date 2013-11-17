#include <iostream>
#include <vector>
using namespace std;

#include "../funopttoolkit/funopttoolkit.hpp"
using namespace funopt;

int main(int argc, char** argv) {
	vector<Triplet> triplets;
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			triplets.push_back(Triplet(i, j, i + j));
		}
	}

	Sparse64f A(triplets, 3, 3);
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			printf("(%d, %d) = %f\n", i, j, A(i, j));
		}
	}

	Sparse64f B = A * A;
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			printf("(%d, %d) = %f\n", i, j, B(i, j));
		}
	}

}
