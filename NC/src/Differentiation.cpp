#include "../include/Differentiation.hpp"
#include <cmath>
#include <stdexcept>

Differentiation::Differentiation(function<double(double)> f, double stepSize)
    : func(f), h(stepSize)
{
    if (!func)
        throw invalid_argument("Function object is empty.");
    if (h <= 0.0)
        throw invalid_argument("Step size h must be positive.");
}
double Differentiation::absoluteError(double exact, double approx) {
    return fabs(exact - approx); // absolute error is always non-negative
}
