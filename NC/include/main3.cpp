#include <iostream>
#include "Root.hpp"

using namespace std;

int main() {
    double a = 0.5, b = 1.0;     // hardcoded interval
    double x0 = 2.0;         
    double dummy = 0.0;

    Root *m1 = new Bisection();
    Root *m2 = new NewtonRaphson();
    Root *m3 = new FixedPoint();

    cout << "Bisection Root      = " << m1->solve(a, b) << endl;
    cout << "Newton Raphson Root = " << m2->solve(x0, dummy) << endl;
    cout << "Fixed Point Root    = " << m3->solve(x0, dummy) << endl;

    delete m1;
    delete m2;
    delete m3;

    return 0;
}









/*
#include <iostream>
#include "root.hpp"
using namespace std;

int main() {
    Root r;

    double a = 0.5, b = 1.0;   // hard-coded interval
    double x0 = 2.0;

    cout << "Bisection Root      : " << r.bisection(a, b) << endl;
    cout << "Newton-Raphson Root : " << r.newtonRaphson(x0) << endl;
    cout << "Fixed Point Root    : " << r.fixedPoint(x0) << endl;

    return 0;
}
    */