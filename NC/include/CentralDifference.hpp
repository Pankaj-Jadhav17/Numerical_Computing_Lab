#ifndef CENTRAL_DIFFERENCE_HPP
#define CENTRAL_DIFFERENCE_HPP
#include "Differentiation.hpp"

// f'(x) ~= [f(x+h) - f(x-h)] / (2h)   -- O(h^2) accurate
class CentralDifference : public Differentiation {
public:
    CentralDifference(function<double(double)> f, double h);
    double derivative(double x) const override;
    string methodName() const override { return "Central Difference"; }
};
#endif
