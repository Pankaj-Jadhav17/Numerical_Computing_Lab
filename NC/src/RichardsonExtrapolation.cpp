#include "../include/RichardsonExtrapolation.hpp"
#include "../include/CentralDifference.hpp"

RichardsonExtrapolation::RichardsonExtrapolation(function<double(double)> f, double h)
    : Differentiation(f, h) {}

double RichardsonExtrapolation::centralAt(double x, double stepSize) const {
    CentralDifference cd(func, stepSize);
    return cd.derivative(x);
}

double RichardsonExtrapolation::derivative(double x) const {
    double D_h    = centralAt(x, h);
    double D_half = centralAt(x, h / 2.0);
    return (4.0 * D_half - D_h) / 3.0;
}
