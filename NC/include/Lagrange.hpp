#ifndef LAGRANGE_HPP
#define LAGRANGE_HPP
#include "Interpolation.hpp"
using namespace std;

class Lagrange : public Interpolation {
public:
    Lagrange(const vector<double>& xs, const vector<double>& ys);
    double interpolate(double x) const override;
    double basisPoly(int i, double x) const;
    string methodName() const override { return "Lagrange Interpolation"; }
};
#endif      
