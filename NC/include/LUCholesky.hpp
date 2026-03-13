#ifndef LU_CHOLESKY_HPP
#define LU_CHOLESKY_HPP

#include "SLE.hpp"

class LUCholesky : public SLE {
public:
    LUCholesky(int r,int c);
    std::vector<double> solve() override;
};

#endif