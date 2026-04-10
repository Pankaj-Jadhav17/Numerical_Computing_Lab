#include "../include/Lagrange.hpp"
#include <stdexcept>
using namespace std;

// Constructor: chain goes  Lagrange --> Interpolation --> Matrix(n, 2)
Lagrange::Lagrange(const vector<double>& xs, const vector<double>& ys) : Interpolation(xs, ys) {}
double Lagrange::basisPoly(int i, double x) const {
    double L = 1.0;
    for (int j = 0; j < rows; ++j) {          // rows = n  (from Matrix)
        if (j == i) continue;
        double xi    = data[i][0];             // data[i][0] = x_i  (from Matrix)
        double xj    = data[j][0];             // data[j][0] = x_j  (from Matrix)
        double denom = xi - xj;
        if (abs(denom) < 1e-15)
            throw runtime_error(
                "Lagrange: duplicate x-values at indices "
                + to_string(i) + " and " + to_string(j));
        L *= (x - xj) / denom;
    }
    return L;
}
double Lagrange::interpolate(double x) const
{
    double result = 0.0;
    for (int i = 0; i < rows; ++i)             // rows = n  (from Matrix)
        result += data[i][1] * basisPoly(i, x); // data[i][1] = y_i  (from Matrix)
    return result;
}