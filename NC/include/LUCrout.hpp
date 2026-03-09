#ifndef LU_CROUT_HPP
#define LU_CROUT_HPP

#include "SLE.hpp"

class LUCrout : public SLE {

public:

    LUCrout(int r,int c);

    std::vector<double> solve() override;
};

#endif