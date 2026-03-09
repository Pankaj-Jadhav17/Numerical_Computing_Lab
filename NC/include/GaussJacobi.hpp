#ifndef GAUSS_JACOBI_HPP
#define GAUSS_JACOBI_HPP

#include "SLE.hpp"

class GaussJacobi : public SLE {

public:
    GaussJacobi(int r, int c);
    std::vector<double> solve() override;
};

#endif