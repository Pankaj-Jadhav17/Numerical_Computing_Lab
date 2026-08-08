#ifndef BACKWARD_DIFFERENCE_HPP
#define BACKWARD_DIFFERENCE_HPP
#include "Differentiation.hpp"

// f'(x) ~= [f(x) - f(x-h)] / h        -- O(h) accurate
class BackwardDifference : public Differentiation {
public:
    BackwardDifference(function<double(double)> f, double h);
    double derivative(double x) const override;
    string methodName() const override { return "Backward Difference"; }
};
#endif
