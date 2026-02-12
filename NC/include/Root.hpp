#ifndef ROOT_HPP
#define ROOT_HPP

class Root {
public:
    Root();  
    virtual ~Root(); 
    virtual double solve(double &x, double &y) = 0;  
};
class Bisection : public Root {  
public: 
    double solve(double &a, double &b) override;  
};
class NewtonRaphson : public Root { 
public:                                            
    double solve(double &x0, double &dummy) override; 
};
class FixedPoint : public Root {
public:
    double solve(double &x0, double &dummy) override;
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