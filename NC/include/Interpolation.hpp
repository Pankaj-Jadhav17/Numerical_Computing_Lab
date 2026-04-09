#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP
#include "Matrix.hpp"
#include <vector>
#include <string>

// Inherits Matrix to store the data point table [x_i | y_i]
class Interpolation : public Matrix {
protected:
    int n;   // number of data points
    std::vector<double> xData;
    std::vector<double> yData;
public:
    // Construct from parallel x / y vectors
    Interpolation(const std::vector<double>& xs, const std::vector<double>& ys);

    virtual ~Interpolation() = default;
    // Evaluate the interpolating polynomial at x
    virtual double interpolate(double x) const = 0;
    // Human-readable method name
    virtual std::string methodName() const = 0;

    // Accessors
    int    numPoints() const { return n; }
    double xAt(int i)  const { return xData[i]; }
    double yAt(int i)  const { return yData[i]; }
};

#endif