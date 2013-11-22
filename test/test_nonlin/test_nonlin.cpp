#include "../../funopttoolkit/funopttoolkit.hpp"
using namespace funopt;

class func : public nonlin::funcNd {
public:
    double operator()(const Vector64f& x) 
    {
        double ret = 0.0;
        for(int i=0; i<x.dim(); i++) {
            double d = x(i) - (i + 1);
            ret += d * d;
        }
        return ret;
    }
};

int main(int argc, char** argv) {
    func f;
    nonlin::Solver solver;
    Vector64f x0 = Vector64f::rand(5);
    Vector64f x_opt;
    solver.solve(&f, x0, x_opt, nonlin::SOLVER_NELDER_SIMPLEX, 200, 1.0e-20);

    cout << x_opt << endl;
    double f_opt = f(x_opt);
    printf("f = %f\n", f_opt);

}
