#ifndef LAGRANGE_HPP
#define LAGRANGE_HPP
#include "Interpolation.hpp"

class Lagrange : public Interpolation {
public:
    Lagrange(const std::vector<double>& xs, const std::vector<double>& ys);
    // Evaluate P(x) using the Lagrange formula
    double interpolate(double x) const override;
    // Evaluate the i-th basis polynomial L_i(x)
    double basisPoly(int i, double x) const;
    std::string methodName() const override { return "Lagrange Interpolation"; }
};
#endif