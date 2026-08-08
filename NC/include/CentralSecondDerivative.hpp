#ifndef CENTRAL_SECOND_DERIVATIVE_HPP
#define CENTRAL_SECOND_DERIVATIVE_HPP
#include "Differentiation.hpp"

// f''(x) ~= [f(x+h) - 2f(x) + f(x-h)] / h^2   -- O(h^2) accurate
// Bonus/extension class: shows the Differentiation hierarchy is not
// limited to first-derivative formulas.
class CentralSecondDerivative : public Differentiation {
public:
    CentralSecondDerivative(function<double(double)> f, double h);
    double derivative(double x) const override;   // returns f''(x), not f'(x)
    string methodName() const override { return "Central Difference (2nd derivative)"; }
};
#endif
