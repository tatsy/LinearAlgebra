#include "../../funopttoolkit/funopttoolkit.hpp"
using namespace funopt;

int main(int argc, char** argv) {
    double elemA[] = {  2.0,  1.0, 1.0, 0.0, 0.0,
                       -1.0, -1.0, 0.0, 1.0, 0.0,
                        1.0,  3.0, 0.0, 0.0, 1.0};
    Matrix64f A(elemA, 3, 5);

    double elemc[] = { -40.0, -60.0, 0.0, 0.0, 0.0 };
    Vector64f c(elemc, 5);

    double elemb[] = { 70.0, 40.0, 90.0 };
    Vector64f b(elemb, 3);

    linprog::Solver solver;
    solver.solve(c, A, b, linprog::SOLVER_PRIMAL_DUAL);
    cout << solver.getXMin() << endl;
    cout << solver.getFMin() << endl;
}
