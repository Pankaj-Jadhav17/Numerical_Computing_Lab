#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP

#include "solver.hpp"

class GaussianElimination : public Solver
{
private:
    bool pivoting;

public:
    GaussianElimination(bool usePivot);
    vector<double> solve(Matrix &A) override;
};

#endif
