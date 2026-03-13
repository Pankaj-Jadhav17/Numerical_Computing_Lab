#ifndef LU_DOOLITTLE_HPP
#define LU_DOOLITTLE_HPP

#include "SLE.hpp"

class LUDoolittle : public SLE {
public:
    LUDoolittle(int r,int c);
    std::vector<double> solve() override;
};

#endif