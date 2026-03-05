#ifndef LU_H
#define LU_H

#include "SLE.hpp"

class LU : public SLE {
public:
    LU(int r, int c);
    void solve(std::ofstream& outFile) override;
};

#endif