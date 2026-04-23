#ifndef LEASTSQUARES_HPP
#define LEASTSQUARES_HPP
#include "CurveFitting.hpp"
#include <vector>
#include <string>
using namespace std;

class LeastSquares : public CurveFitting {
public:
    explicit LeastSquares(const Matrix& m);  // delegates validation to CurveFitting
    //  virtual ~LeastSquares() = default;         // default destructor
    // Solves the normal equations and returns {a, b}  (a = slope, b = intercept).
    vector<double> fit() const override;

    double evaluate(double x, const vector<double>& coeffs) const override;
    double rmsError(const vector<double>& coeffs) const override;

    string methodName() const override { return "Least Squares Linear Fit"; }

    double evaluate(double x, double a, double b) const;
    double rmsError(double a, double b) const;
};
#endif
