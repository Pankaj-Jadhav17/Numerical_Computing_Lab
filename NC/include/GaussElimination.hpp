#ifndef GAUSS_ELIMINATION_HPP
#define GAUSS_ELIMINATION_HPP

#include "SLE.hpp"

class GaussElimination : public SLE {

public:

    GaussElimination(int r,int c);

    std::vector<double> solve() override;
};

#endif 