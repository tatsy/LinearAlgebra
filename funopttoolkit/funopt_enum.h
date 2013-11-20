#ifndef _FUNOPT_ENUM_H_
#define _FUNOPT_ENUM_H_

namespace funopt {
    enum SolverType {
        SOLVER_LU = 0x01,
	    SOLVER_QR = 0x02,
        SOLVER_CG = 0x04
    };
}

#endif
