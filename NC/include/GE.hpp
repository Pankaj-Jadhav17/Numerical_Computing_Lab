#ifndef GE_H
#define GE_H

#include "SLE.hpp"

class GE : public SLE
{
private:
    bool pivoting;

public:
    GE(int n, bool usePivot);

    vector<double> solve() override;
};

#endif