#ifndef ROOT_HPP
#define ROOT_HPP

class Root {
public:
    Root();  // constructor
    virtual ~Root(); // virtual destructor to ensure proper cleanup of derived class objects through base class pointers

    virtual double solve(double &x, double &y) = 0;  // pure virtual function to be overridden by derived classes
};

class Bisection : public Root {  // inherits from Root
public:  // override the solve method for bisection method
    double solve(double &a, double &b) override;  // override keyword indicates this method overrides a virtual method in the base class
};

// Newton-Raphson method class
class NewtonRaphson : public Root { 
public:                                             // override the solve method for Newton-Raphson method
    double solve(double &x0, double &dummy) override; // dummy parameter is not used but required to match the base class method signature
};

class FixedPoint : public Root {
public:
    double solve(double &x0, double &dummy) override;// override the solve method for fixed-point iteration method
};

#endif


/*
#ifndef ROOT_HPP
#define ROOT_HPP

class Root {
public:
    Root();
    ~Root();

    double bisection(double &a, double &b);
    double newtonRaphson(double &x0);
    double fixedPoint(double &x0);
};

#endif
*/