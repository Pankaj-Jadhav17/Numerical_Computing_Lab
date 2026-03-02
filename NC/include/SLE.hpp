#ifndef SLE_H
#define SLE_H

#include "Matrix.hpp"

class SLE : public Matrix
{
public:
    SLE(int n = 0);

    // Pure virtual solver
    virtual vector<double> solve() = 0;

    virtual ~SLE() {}
};

#endif