#ifndef FORWARD_DIFFERENCE_HPP
#define FORWARD_DIFFERENCE_HPP
#include "Differentiation.hpp"

// f'(x) ~= [f(x+h) - f(x)] / h        -- O(h) accurate
class ForwardDifference : public Differentiation {
public:
    ForwardDifference(function<double(double)> f, double h);
    double derivative(double x) const override;
    string methodName() const override { return "Forward Difference"; }
};
#endif
