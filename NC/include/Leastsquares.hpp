#ifndef LEASTSQUARES_HPP
#define LEASTSQUARES_HPP

#include "Matrix.hpp"
#include <string>
using namespace std;

// LeastSquares is CHILD of Matrix
//
// Fits a straight line:   y = a*x + b
//
// Normal equations (from minimizing S = Σ(fi - (a*xi + b))²):
//
//   [ Σxi²   Σxi  ] [a]   [Σfi*xi]
//   [ Σxi    N    ] [b] = [Σfi   ]
//
// Matrix data layout (n rows, 2 cols):
//   data[i][0] = xi   (x value of i-th point)
//   data[i][1] = fi   (y/f value of i-th point)
//   rows       = N    (number of data points)
//   cols       = 2    (always)
//
// NO own data members — uses only Matrix members: rows, cols, data

class LeastSquares : public Matrix {

public:
    // Constructor: stores (xi, fi) pairs into Matrix data[][]
    LeastSquares(const vector<double>& xs, const vector<double>& fs);

    // Solves the 2x2 normal equation system
    // Returns {a, b} where y = a*x + b
    // Uses Matrix::inverse() to solve [A]{a,b} = {rhs}
    vector<double> fit() const;

    // Evaluates fitted line at a given x:  y = a*x + b
    double evaluate(double x, double a, double b) const;

    // Computes RMS error = sqrt( (1/N) * Σ(fi - (a*xi+b))² )
    // From notes: rms = sqrt( S / N )  where S = Σ(yi-fi)²
    double rmsError(double a, double b) const;

    // Accessors — read from Matrix members only
    int    numPoints() const { return rows;       }   // rows = N
    double xAt(int i)  const { return data[i][0]; }   // col 0 = xi
    double fAt(int i)  const { return data[i][1]; }   // col 1 = fi

    string methodName() const { return "Least Squares Linear Fit"; }
};

#endif