#include "../include/CentralSecondDerivative.hpp"

CentralSecondDerivative::CentralSecondDerivative(function<double(double)> f, double h)
    : Differentiation(f, h) {}

double CentralSecondDerivative::derivative(double x) const {
    return (func(x + h) - 2.0 * func(x) + func(x - h)) / (h * h);
}
