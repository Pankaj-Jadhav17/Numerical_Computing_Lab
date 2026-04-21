#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "Matrix.hpp"
#include <string>
using namespace std;

class Interpolation : public Matrix {

public:
    // Constructor: takes parallel x and y vectors
    Interpolation(const vector<double>& xs, const vector<double>& ys);

    virtual ~Interpolation() = default;

    // Pure virtual — subclass MUST implement
    virtual double interpolate(double x) const = 0;
    virtual string methodName()          const = 0;

    // Accessors — all use Matrix members (rows, data)
    int    numPoints() const { return rows;         }  // rows  = n
    double xAt(int i)  const { return data[i][0];   }  // col 0 = x
    double yAt(int i)  const { return data[i][1];   }  // col 1 = y
};

#endif