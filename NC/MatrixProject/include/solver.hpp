#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "matrix.hpp"
#include <vector>

class Solver
{
public:
    virtual vector<double> solve(Matrix &A) = 0;
    virtual ~Solver() {}
};

#endif
