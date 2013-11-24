#include "../../funopttoolkit/funopttoolkit.hpp"
#include "../../funopttoolkit/Brent.h"
using namespace funopt;

double f(const Vector64f& x) 
{
    double ret = 0.0;
    for(int i=0; i<x.dim(); i++) {
        double d = x(i) - (i + 1);
        ret += d * d;
    }
    return ret;
}

double g(double x)
{
    return (x - 1.0) * (x - 3.0);
}

int main(int argc, char** argv) {
    nonlin::funcNd func(f);

    nonlin::Solver solver;
    Vector64f x0 = Vector64f::rand(5);
    Vector64f x_opt;
    solver.solve(func, x0, x_opt, nonlin::SOLVER_POWELL, 200, 1.0e-20);

    cout << x_opt << endl;
    double f_opt = func(x_opt);
    printf("f = %f\n", f_opt);
}
