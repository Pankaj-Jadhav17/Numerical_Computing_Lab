#ifndef GAUSS_ELIMINATION_H
#define GAUSS_ELIMINATION_H

#include "SLE.hpp"

class GaussElimination : public SLE {
private:
    bool usePivoting; 

public:
    GaussElimination(int r, int c, bool pivot);
    void solve(std::ofstream& outFile) override;
};

#endif 