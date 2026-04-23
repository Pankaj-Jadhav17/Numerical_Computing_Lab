#include "../include/LeastSquares.hpp"
#include <stdexcept>
#include <cmath>
using namespace std;

LeastSquares::LeastSquares(const Matrix& m) : CurveFitting(m) {}
vector<double> LeastSquares::fit() const
{
    double sumX  = 0.0;
    double sumX2 = 0.0;
    double sumF  = 0.0;
    double sumFX = 0.0;
    for (int i = 0; i < rows; ++i) {
        double xi = data[i][0];
        double fi = data[i][1];
        sumX  += xi;
        sumX2 += xi * xi;
        sumF  += fi;
        sumFX += fi * xi;
    }

    double det = sumX2 * static_cast<double>(rows) - sumX * sumX;
    if (abs(det) < 1e-12)
        throw runtime_error(
            "LeastSquares: System is singular. "
            "Check for duplicate or collinear x-values.");

    double a = (sumFX * static_cast<double>(rows) - sumF * sumX) / det;
    double b = (sumX2 * sumF - sumX * sumFX) / det;

    return {a, b};
}
double LeastSquares::evaluate(double x, const vector<double>& coeffs) const {
    if (coeffs.size() < 2)
        throw invalid_argument(
            "LeastSquares::evaluate: coeffs must contain at least {a, b}.");
    return coeffs[0] * x + coeffs[1];
}

double LeastSquares::rmsError(const vector<double>& coeffs) const {
    if (coeffs.size() < 2)
        throw invalid_argument(
            "LeastSquares::rmsError: coeffs must contain at least {a, b}.");
    if (rows < 1)
        return 0.0;

    double S = 0.0;
    for (int i = 0; i < rows; ++i) {
        double xi      = data[i][0];
        double fi      = data[i][1];
        double residual = fi - (coeffs[0] * xi + coeffs[1]);
        S += residual * residual;
    }
    return sqrt(S / static_cast<double>(rows));
}

double LeastSquares::evaluate(double x, double a, double b) const {
    return evaluate(x, vector<double>{a, b});
}

double LeastSquares::rmsError(double a, double b) const {
    return rmsError(vector<double>{a, b});
}