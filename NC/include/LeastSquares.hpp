#ifndef LEASTSQUARES_HPP
#define LEASTSQUARES_HPP

#include <vector>
#include <string>
using namespace std;

// LeastSquares: Fits a straight line y = a*x + b to data points
// Using least squares minimization: minimize S = Σ(fi - (a*xi + b))²
//
// Normal equations derived from ∂S/∂a = 0 and ∂S/∂b = 0:
//   [ Σxi²   Σxi  ] [a]   [Σfi*xi]
//   [ Σxi    N    ] [b] = [Σfi   ]

class LeastSquares {
private:
    vector<double> xs;   // x values
    vector<double> fs;   // y/f values
    int n;               // number of data points

public:
    // Constructor: stores the data points
    LeastSquares(const vector<double>& xs_in, const vector<double>& fs_in);

    // Solves the 2x2 normal equation system
    // Returns {a, b} where y = a*x + b
    vector<double> fit() const;

    // Evaluates fitted line at a given x: y = a*x + b
    double evaluate(double x, double a, double b) const;

    // Computes RMS error = sqrt( (1/N) * Σ(fi - (a*xi+b))² )
    double rmsError(double a, double b) const;

    // Accessors
    int    numPoints() const { return n; }
    double xAt(int i) const  { return xs[i]; }
    double fAt(int i) const  { return fs[i]; }

    string methodName() const { return "Least Squares Linear Fit"; }
};

#endif
