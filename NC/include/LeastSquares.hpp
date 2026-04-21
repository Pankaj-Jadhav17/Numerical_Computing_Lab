#ifndef LEASTSQUARES_HPP
#define LEASTSQUARES_HPP
#include "Matrix.hpp"      // inherits rows, cols, data
#include <vector>
#include <string>
using namespace std;

class LeastSquares : public Matrix {

public:
    explicit LeastSquares(const Matrix& m);
    vector<double> fit() const;

    // Evaluates fitted line at a given x: returns a*x + b
    double evaluate(double x, double a, double b) const;

    // Computes RMS error = sqrt( (1/N) * Σ(fi - (a*xi+b))² )
    double rmsError(double a, double b) const;

    //Accessors (use Matrix's protected members, no new storage) 
    int    numPoints() const { return rows; }          // rows  == n
    double xAt(int i)  const { return data[i][0]; }   // col 0 == xi
    double fAt(int i)  const { return data[i][1]; }   // col 1 == fi

    string methodName() const { return "Least Squares Linear Fit"; }
};
#endif