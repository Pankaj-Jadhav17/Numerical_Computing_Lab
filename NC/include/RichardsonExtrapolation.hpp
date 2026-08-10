#ifndef RICHARDSON_EXTRAPOLATION_HPP
#define RICHARDSON_EXTRAPOLATION_HPP
#include "Differentiation.hpp"

class RichardsonExtrapolation : public Differentiation {
public:
    RichardsonExtrapolation(function<double(double)> f, double h);
    double derivative(double x) const override;

    // wants both D(h) (Central Difference) and R(h) reported side by side.
    double centralAt(double x, double stepSize) const;

    string methodName() const override { return "Richardson Extrapolation"; }
};
#endif
