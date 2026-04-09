#include "../include/Lagrange.hpp"
#include <stdexcept>
using namespace std;          

Lagrange::Lagrange(const vector<double>& xs,const vector<double>& ys): Interpolation(xs, ys) {}

double Lagrange::basisPoly(int i, double x) const
{
    double L = 1.0;
    for (int j = 0; j < n; ++j) {
        if (j == i) continue;
        double denom = xData[i] - xData[j];
        if (abs(denom) < 1e-15)              
            throw runtime_error(
                "Lagrange: duplicate x-values at indices "
                + to_string(i) + " and " + to_string(j));
        L *= (x - xData[j]) / denom;
    }
    return L;
}
double Lagrange::interpolate(double x) const
{
    double result = 0.0;
    for (int i = 0; i < n; ++i)
        result += yData[i] * basisPoly(i, x);
    return result;
}