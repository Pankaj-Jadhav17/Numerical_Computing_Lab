#include "Root.hpp"
#include <cmath>

// constructor & destructor
Root::Root() {}
Root::~Root() {}

// Example function to find root of: f(x) = 4x³ - 3x
static double f(double x) {
    return 4*x*x * x - 3*x;
}

static double df(double x) {
    return 12 * x * x - 3;
}

// g(x) for fixed point
static double g(double x) {
    return sqrt(2);  
}

//Bisection 
double Bisection::solve(double &a, double &b) {
    double m;
    for (int i = 0; i < 20; i++) {
        m = (a + b) / 2;
        if (f(a) * f(m) < 0)
            b = m;
        else
            a = m;
    }
    return m;
}

//Newton Raphson 
double NewtonRaphson::solve(double &x0, double &dummy) {
    for (int i = 0; i < 10; i++)
        x0 = x0 - f(x0) / df(x0);     
    return x0;
}

// Fixed Point Iteration
double FixedPoint::solve(double &x0, double &dummy) {
    for (int i = 0; i < 10; i++)
        x0 = g(x0);
    return x0;
}





























/*
#include "root.hpp"
#include <cmath>

// Constructor & Destructor
Root::Root() {}
Root::~Root() {}

// f(x) = 4x³ - 3x (example function to find root of) - can be modified as needed
static double f(double x) {
    return 4 * x * x * x - 3 * x;
}

// f'(x) = 12x² - 3
static double df(double x) {
    return 12 * x * x - 3;
}

// g(x) for fixed point iteration
static double g(double x) {
    return sqrt(2);
}
// Bisection Method
double Root::bisection(double &a, double &b) {
    double m;
    for (int i = 0; i < 20; i++) {
        m = (a + b) / 2;
        if (f(a) * f(m) < 0)
            b = m;
        else
            a = m;
    }
    return m;
}
// Newton–Raphson Method
double Root::newtonRaphson(double &x0) {
    for (int i = 0; i < 10; i++)
        x0 = x0 - f(x0) / df(x0);
    return x0;
}
// Fixed Point Iteration
double Root::fixedPoint(double &x0) {
    for (int i = 0; i < 10; i++)
        x0 = g(x0);
    return x0;
}

 */