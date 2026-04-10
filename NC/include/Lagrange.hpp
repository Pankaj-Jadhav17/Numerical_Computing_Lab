#ifndef LAGRANGE_HPP
#define LAGRANGE_HPP
#include "Interpolation.hpp"
using namespace std;

class Lagrange : public Interpolation {

public:
    // Constructor: passes xs, ys up to Interpolation --> Matrix
    Lagrange(const vector<double>& xs, const vector<double>& ys);
    // Evaluates P(x) = sum of y_i * L_i(x)
    double interpolate(double x) const override;
    // Evaluates i-th Lagrange basis polynomial:
    double basisPoly(int i, double x) const;
    // Implements pure virtual from Interpolation
    string methodName() const override { return "Lagrange Interpolation"; }
};

#endif      