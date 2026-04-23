#ifndef CURVEFITTING_HPP
#define CURVEFITTING_HPP
#include "Matrix.hpp"
#include <vector>
#include <string>
using namespace std;

class CurveFitting : public Matrix {
public:
    explicit CurveFitting(const Matrix& m);   // defined in CurveFitting.cpp
    virtual ~CurveFitting();                  // defined in CurveFitting.cpp

    //Pure virtual interface
    virtual vector<double> fit() const = 0;
    virtual double evaluate(double x, const vector<double>& coeffs) const = 0;
    virtual double rmsError(const vector<double>& coeffs) const = 0;
    virtual string methodName() const = 0;

    // All subclasses read the same 2-column data layout (col-0 = x, col-1 = f).
    int    numPoints() const { return rows; }
    double xAt(int i)  const { return data[i][0]; }
    double fAt(int i)  const { return data[i][1]; }
};
#endif 