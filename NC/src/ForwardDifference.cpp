#include "../include/ForwardDifference.hpp"

ForwardDifference::ForwardDifference(function<double(double)> f, double h)
    : Differentiation(f, h) {}

double ForwardDifference::derivative(double x) const {
    return (func(x + h) - func(x)) / h;
}
