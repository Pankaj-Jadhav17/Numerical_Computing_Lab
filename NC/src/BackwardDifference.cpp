#include "../include/BackwardDifference.hpp"

BackwardDifference::BackwardDifference(function<double(double)> f, double h)
    : Differentiation(f, h) {}

double BackwardDifference::derivative(double x) const {
    return (func(x) - func(x - h)) / h;
}
