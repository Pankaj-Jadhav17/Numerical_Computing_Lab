#include "../include/CentralDifference.hpp"

CentralDifference::CentralDifference(function<double(double)> f, double h)
    : Differentiation(f, h) {}

double CentralDifference::derivative(double x) const {
    return (func(x + h) - func(x - h)) / (2.0 * h);
}
