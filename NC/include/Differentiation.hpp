#ifndef DIFFERENTIATION_HPP
#define DIFFERENTIATION_HPP
#include <functional>
#include <string>
using namespace std;

class Differentiation {
protected:
    function<double(double)> func;   // the function f(x) being differentiated
    double h;                        // step size

public:
    Differentiation(function<double(double)> f, double stepSize);
    virtual ~Differentiation() = default;

    // Pure virtual — subclass MUST implement its own formula
    virtual double derivative(double x) const = 0;
    virtual string methodName() const = 0;

    // Accessors
    double getStepSize() const { return h; }
    void   setStepSize(double newH) { h = newH; }

    // Shared utility used by every subclass / by the analyzer
    static double absoluteError(double exact, double approx);
};

#endif
