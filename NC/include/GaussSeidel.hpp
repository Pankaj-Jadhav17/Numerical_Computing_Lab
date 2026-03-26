#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP

#include "SLE.hpp"

class GaussSeidel : public SLE {

public:
    GaussSeidel(int r, int c);
    std::vector<double> solve() override;
};

#endif