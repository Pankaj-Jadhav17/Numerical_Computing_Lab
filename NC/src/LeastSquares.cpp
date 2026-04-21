#include "../include/LeastSquares.hpp"
#include <stdexcept>
#include <cmath>
using namespace std;

LeastSquares::LeastSquares(const Matrix& m) : Matrix(m) { 
    if (cols != 2)
        throw invalid_argument(
            "LeastSquares: Matrix must have exactly 2 columns "
            "(col 0 = x values, col 1 = f values).");
    if (rows < 2)
        throw invalid_argument(
            "LeastSquares: Need at least 2 data points (rows).");
}

vector<double> LeastSquares::fit() const
{
    // Accumulate sums using inherited data
    double sumX  = 0.0;   // Σxi
    double sumX2 = 0.0;   // Σxi²
    double sumF  = 0.0;   // Σfi
    double sumFX = 0.0;   // Σfi·xi

    for (int i = 0; i < rows; ++i) {       // rows == N (inherited)
        double xi = data[i][0];             // inherited data col-0 = xi
        double fi = data[i][1];             // inherited data col-1 = fi
        sumX  += xi;
        sumX2 += xi * xi;
        sumF  += fi;
        sumFX += fi * xi;
    }

    double det = sumX2 * (double)rows - sumX * sumX;

    if (abs(det) < 1e-12)
        throw runtime_error(
            "LeastSquares: System is singular. "
            "Check for duplicate or collinear x-values.");

    // Cramer's rule
    double a = (sumFX * (double)rows - sumF * sumX) / det;
    double b = (sumX2 * sumF          - sumX * sumFX) / det;

    return {a, b};
}

double LeastSquares::evaluate(double x, double a, double b) const
{
    return a * x + b;
}

double LeastSquares::rmsError(double a, double b) const {
    if (rows < 1)
        return 0.0;

    double S = 0.0;

    for (int i = 0; i < rows; ++i) {           // rows == N (inherited)
        double xi      = data[i][0];            // xi from inherited data
        double fi      = data[i][1];            // fi from inherited data
        double fi_fit  = a * xi + b;            // fitted value
        double residual = fi - fi_fit;
        S += residual * residual;
    }

    return sqrt(S / (double)rows);
}