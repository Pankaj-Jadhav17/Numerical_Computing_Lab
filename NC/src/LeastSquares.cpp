#include "../include/LeastSquares.hpp"
#include <stdexcept>
#include <cmath>
using namespace std;

// Constructor: stores the data points
LeastSquares::LeastSquares(const vector<double>& xs_in, const vector<double>& fs_in)
    : xs(xs_in), fs(fs_in), n((int)xs_in.size())
{
    if (xs.size() != fs.size() || xs.empty())
        throw invalid_argument("LeastSquares: x and f vectors must be non-empty and same length.");
}

// fit(): Builds and solves the 2x2 normal equation system
//
// From least squares minimization (minimize S = Σ(fi - (a*xi + b))²):
//   ∂S/∂a = 0  =>  a·Σxi²  + b·Σxi  = Σfi·xi
//   ∂S/∂b = 0  =>  a·Σxi   + b·N    = Σfi
//
// Matrix form [A]{coeff} = {rhs}:
//   A    = [ Σxi²  Σxi ]
//          [ Σxi   N   ]
//   rhs  = [ Σfi·xi ]
//          [ Σfi    ]
//
// Uses Cramer's rule to solve the 2x2 system

vector<double> LeastSquares::fit() const
{
    if (n < 2)
        throw invalid_argument("LeastSquares: Need at least 2 data points.");

    // Accumulate sums
    double sumX  = 0.0;   // Σxi
    double sumX2 = 0.0;   // Σxi²
    double sumF  = 0.0;   // Σfi
    double sumFX = 0.0;   // Σfi·xi

    for (int i = 0; i < n; ++i) {
        double xi = xs[i];
        double fi = fs[i];
        sumX  += xi;
        sumX2 += xi * xi;
        sumF  += fi;
        sumFX += fi * xi;
    }

    // Build the 2x2 system using Cramer's rule
    // A = [ Σxi²  Σxi ]    rhs = [ Σfi·xi ]
    //     [ Σxi   N   ]          [ Σfi    ]

    double det = sumX2 * (double)n - sumX * sumX;

    if (abs(det) < 1e-12)
        throw runtime_error("LeastSquares: System is singular. "
                            "Check for duplicate or collinear x-values.");

    // Solve using Cramer's rule
    double a = (sumFX * (double)n - sumF * sumX) / det;
    double b = (sumX2 * sumF - sumX * sumFX) / det;

    return {a, b};
}

// evaluate(): Returns y = a*x + b at given x
double LeastSquares::evaluate(double x, double a, double b) const
{
    return a * x + b;
}

// rmsError(): Computes Root Mean Square error
// RMS = sqrt( (1/N) * Σ(fi - fitted_yi)² )

double LeastSquares::rmsError(double a, double b) const
{
    if (n < 1)
        return 0.0;

    double S = 0.0;

    for (int i = 0; i < n; ++i) {
        double xi      = xs[i];
        double yi      = fs[i];              // actual y value
        double yi_fit  = a * xi + b;         // fitted value
        double residual = yi - yi_fit;
        S += residual * residual;            // S = Σ(yi - fi)²
    }

    // RMS = sqrt( S / N )
    return sqrt(S / (double)n);
}
