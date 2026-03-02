#include "BM.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;
double BisectionMethod::f(double x) {
    return 4*x*x*x - 3*x;
}
BisectionMethod::BisectionMethod(double left, double right, int tolPower)
    : a(left), b(right), k(tolPower) {
    tolerance = pow(10, -k);
}

bool BisectionMethod::validateInterval() {
    return f(a) * f(b) < 0;
}

double BisectionMethod::findRoot() {
    history.clear();

    double left = a, right = b;
    double mid, prev = left;

    for (int i = 1; i <= 1000; i++) {
        mid = (left + right) / 2.0;
        double err = fabs(mid - prev);

        history.push_back({
            i, left, right, mid, f(mid), err, right - left
        });
        if (i > 1 && err < tolerance)
            return mid;
        if (f(left) * f(mid) < 0)
            right = mid;
        else
            left = mid;
        prev = mid;
    }
    return mid;
}
void BisectionMethod::displayIterationTable() {
    cout << "\nIter      a           b           mid        f(mid)        error\n";
    cout << "------------------------------------------------------------------\n";
    for (auto &d : history) {
        cout << setw(3) << d.iteration
             << setw(12) << fixed << setprecision(6) << d.a
             << setw(12) << d.b
             << setw(12) << d.midpoint
             << setw(14) << scientific << d.f_midpoint
             << setw(12) << d.error << endl;
    }
}
double BisectionMethod::getTolerance() const {
    return tolerance;
}
int BisectionMethod::getActualIterations() const {
    return history.size();
}


/*
    BM.cpp - Implementation File for Bisection Method
    Function: f(x) = 4x³ - 3x
*/
// #include "BM.hpp"
// #include <iostream>
// #include <cmath>
// #include <iomanip>

// using namespace std;

// // Define the function f(x) = 4x³ - 3x
// double BisectionMethod::f(double x) {
//     return 4 * x * x * x - 3 * x;
// }

// // Constructor
// BisectionMethod::BisectionMethod(double left, double right, int tolerance_power)
//     : a(left), b(right), k(tolerance_power) {
//     tolerance = pow(10, -k);
//     history.clear();
// }

// // Validate that the interval brackets a root
// bool BisectionMethod::validateInterval() {
//     double fa = f(a);
//     double fb = f(b);
    
//     cout << "\n Validating interval [" << a << ", " << b << "]\n";
//     cout << "f(" << a << ") = " << fa << endl;
//     cout << "f(" << b << ") = " << fb << endl;
    
//     if (fa * fb >= 0) {
//         cerr << "\n ERROR: f(a) and f(b) must have opposite signs!" << endl;
//         cerr << "The interval [" << a << ", " << b << "] does not bracket a root.\n";
//         return false;
//     }
    
//     cout << " Valid interval - root exists in [" << a << ", " << b << "]\n";
//     return true;
// }

// // Step 2: Calculate theoretical number of iterations
// int BisectionMethod::calculateTheoreticalIterations() {
//     double n = (log(b - a) - log(tolerance)) / log(2);
//     return (int)ceil(n);
// }

// // Step 1: Find the root using bisection method
// double BisectionMethod::findRoot() {
//     // Clear previous history
//     history.clear();
    
//     double m_prev = a;      // Previous midpoint
//     double m_curr;          // Current midpoint
//     double left = a;        // Working left endpoint
//     double right = b;       // Working right endpoint
//     int iterations = 0;
    
//     cout<<"Bisection method iterations:\n";
//     while (iterations < 1000) {
//         // Calculate midpoint
//         m_curr = (left + right) / 2.0;
//         iterations++;
        
//         // Calculate function value at midpoint
//         double f_mid = f(m_curr);
        
//         // Calculate error
//         double error = fabs(m_curr - m_prev);
        
//         // Calculate interval length
//         double interval_len = right - left;
        
//         // Store iteration data
//         IterationData data;
//         data.iteration = iterations;
//         data.a = left;
//         data.b = right;
//         data.midpoint = m_curr;
//         data.f_midpoint = f_mid;
//         data.error = error;
//         data.interval_length = interval_len;
//         history.push_back(data);
        
//         // Check stopping condition
//         if (error < tolerance && iterations > 1) {
//             return m_curr;
//         }
        
//         // Update interval based on sign test
//         double f_left = f(left);
        
//         if (f_left * f_mid < 0) {
//             // Root is in left half
//             right = m_curr;
//         } else {
//             // Root is in right half
//             left = m_curr;
//         }
        
//         // Update previous midpoint
//         m_prev = m_curr;
//     }
    
//     cerr << "Warning: Maximum iterations reached!\n";
//     return m_curr;
// }

// // Step 4: Get iteration history
// const vector<IterationData>& BisectionMethod::getIterationHistory() const {
//     return history;
// }

// // Step 5: Display iteration table showing interval reduction
// void BisectionMethod::displayIterationTable() {
//     cout << "\n Iteration Table\n";
//     cout << fixed;
//     cout << setw(6) << "Iter" 
//          << setw(15) << "a (Left)" 
//          << setw(15) << "b (Right)"
//          << setw(16) << "Midpoint (m)"
//          << setw(15) << "f(m)"
//          << setw(15) << "Error"
//          << setw(14) << "Interval" << endl;
//     cout << string(96, '-') << endl;
    
//     for (const auto& data : history) {
//         cout << setw(6) << data.iteration
//              << setw(15) << setprecision(8) << data.a
//              << setw(15) << setprecision(8) << data.b
//              << setw(16) << setprecision(10) << data.midpoint
//              << setw(15) << scientific << setprecision(4) << data.f_midpoint
//              << setw(15) << scientific << setprecision(4) << data.error
//              << setw(14) << fixed << setprecision(6) << data.interval_length << endl;
//     }
    
//     cout << string(96, '-') << endl;
// }

// // Step 6: Display final results including tolerance
// void BisectionMethod::displayResults(double root, int actualIterations, int theoreticalIterations) {
//     cout << "\n Final Results\n";
    
//     cout << fixed << setprecision(12);
//     cout << "\n1. ROOT FOUND:\n";
//     cout << "   Approximate root (x):     " << root << endl;
//     cout << "   Function value f(x):      " << scientific << setprecision(6) << f(root) << endl;
    
//     if (!history.empty()) {
//         auto& lastIter = history.back();
//         cout << fixed << setprecision(10);
//         cout << "   Final interval:           [" << lastIter.a << ", " << lastIter.b << "]\n";
//         cout << "   Final interval length:    " << scientific << setprecision(6) << lastIter.interval_length << endl;
//     }
    
//     cout << "\n2. NUMBER OF ITERATIONS:\n";
//     cout << "   Actual iterations:        " << actualIterations << endl;
//     cout << "   Theoretical iterations:   " << theoreticalIterations << endl;
//     cout << "   Difference:               " << (actualIterations - theoreticalIterations) << endl;
    
//     cout << "\n3. TOLERANCE:\n";
//     cout << "   Required tolerance:       " << scientific << setprecision(6) << tolerance << endl;
//     cout << "   Tolerance power (k):      " << k << " (10^-" << k << ")" << endl;
//     if (!history.empty()) {
//         cout << "   Achieved error:           " << scientific << setprecision(6) << history.back().error << endl;
//     }
    
//     cout << "\n4. VERIFICATION:\n";
//     cout << "   Function: f(x) = 4x³ - 3x\n";
//     cout << "   Factored: f(x) = x(4x² - 3)\n";
//     cout << "   Known roots: x = 0, x = ±√(3/4) ≈ ±0.866025\n";
    
//     // Identify which root was found
//     if (fabs(root) < 0.1) {
//         cout << "   ✓ Root found: x = 0\n";
//     } else if (root > 0) {
//         double exact = sqrt(3.0 / 4.0);
//         cout << "   ✓ Root found: x = √(3/4) (positive root)\n";
//         cout << "   Exact value:              " << fixed << setprecision(12) << exact << endl;
//         cout << "   Absolute error:           " << scientific << setprecision(6) << fabs(root - exact) << endl;
//     } else {
//         double exact = -sqrt(3.0 / 4.0);
//         cout << "   Root found: x = -√(3/4) (negative root)\n";
//         cout << "   Exact value:              " << fixed << setprecision(12) << exact << endl;
//         cout << "   Absolute error:           " << scientific << setprecision(6) << fabs(root - exact) << endl;
//     }
//     cout << "\n";
// }

// // Get tolerance value
// double BisectionMethod::getTolerance() const {
//     return tolerance;
// }

// // Get actual number of iterations
// int BisectionMethod::getActualIterations() const {
//     return history.size();
// }

// // Evaluate the function f(x) = 4x³ - 3x
// double BisectionMethod::evaluateFunction(double x) const {
//     return 4 * x * x * x - 3 * x;
// }



























// #include "BM.hpp"
// #include <iostream>
// #include <cmath>
// #include <iomanip>

// using namespace std;

// // Constructor
// BisectionMethod::BisectionMethod(std::function<double(double)> func, 
//                                  double left, 
//                                  double right, 
//                                  double tol,
//                                  int maxIter)
//     : function(func), a(left), b(right), tolerance(tol), maxIterations(maxIter) {
//     history.clear();
// }

// // Validate that the interval brackets a root
// bool BisectionMethod::validateInterval() const {
//     double fa = function(a);
//     double fb = function(b);
    
//     if (fa * fb >= 0) {
//         cerr << "Error: f(a) and f(b) must have opposite signs!" << endl;
//         cerr << "f(" << a << ") = " << fa << endl;
//         cerr << "f(" << b << ") = " << fb << endl;
//         return false;
//     }
//     return true;
// }

// // Calculate theoretical minimum number of iterations
// int BisectionMethod::getTheoreticalIterations() const {
//     return (int)ceil((log(b - a) - log(tolerance)) / log(2));
// }

// // Main bisection algorithm
// double BisectionMethod::findRoot() {
//     // Validate interval first
//     if (!validateInterval()) {
//         return NAN;  // Return Not-a-Number if invalid
//     }
    
//     // Clear previous history
//     history.clear();
    
//     double m_prev = a;  // Previous midpoint
//     double m_curr;      // Current midpoint
//     double left = a;    // Working left endpoint
//     double right = b;   // Working right endpoint
//     int iterations = 0;
    
//     while (iterations < maxIterations) {
//         // Calculate midpoint
//         m_curr = (left + right) / 2.0;
//         iterations++;
        
//         // Calculate error
//         double error = fabs(m_curr - m_prev);
        
//         // Store iteration data
//         IterationData data;
//         data.iteration = iterations;
//         data.midpoint = m_curr;
//         data.error = error;
//         data.intervalLeft = left;
//         data.intervalRight = right;
//         history.push_back(data);
        
//         // Check stopping condition
//         if (error < tolerance && iterations > 1) {
//             return m_curr;
//         }
        
//         // Update interval based on sign test
//         double f_left = function(left);
//         double f_mid = function(m_curr);
        
//         if (f_left * f_mid < 0) {
//             // Root is in left half
//             right = m_curr;
//         } else {
//             // Root is in right half
//             left = m_curr;
//         }
        
//         // Update previous midpoint
//         m_prev = m_curr;
//     }
    
//     // Maximum iterations reached
//     cerr << "Warning: Maximum iterations (" << maxIterations << ") reached!" << endl;
//     return m_curr;
// }

// // Display iteration table
// void BisectionMethod::displayIterationTable() const {
//     cout << "\n" << string(75, '=') << endl;
//     cout << "ITERATION DETAILS" << endl;
//     cout << string(75, '=') << endl;
    
//     cout << setw(6) << "Iter" 
//          << setw(15) << "Midpoint" 
//          << setw(15) << "Error"
//          << setw(15) << "Left [a]"
//          << setw(15) << "Right [b]" << endl;
//     cout << string(75, '-') << endl;
    
//     for (const auto& data : history) {
//         cout << setw(6) << data.iteration
//              << setw(15) << fixed << setprecision(10) << data.midpoint
//              << setw(15) << scientific << setprecision(4) << data.error
//              << setw(15) << fixed << setprecision(6) << data.intervalLeft
//              << setw(15) << fixed << setprecision(6) << data.intervalRight
//              << endl;
//     }
//     cout << string(75, '=') << endl;
// }

// // Display final results
// void BisectionMethod::displayResults(double root, int actualIterations) const {
//     cout << "\n" << string(50, '=') << endl;
//     cout << "BISECTION METHOD RESULTS" << endl;
//     cout << string(50, '=') << endl;
//     cout << fixed << setprecision(12);
//     cout << "Approximate root:        " << root << endl;
//     cout << "Function value f(root):  " << scientific << setprecision(6) 
//          << function(root) << endl;
//     cout << fixed << setprecision(10);
//     cout << "Tolerance:               " << tolerance << endl;
//     cout << "Actual iterations:       " << actualIterations << endl;
//     cout << "Theoretical iterations:  " << getTheoreticalIterations() << endl;
//     cout << string(50, '=') << endl;
// }

// // Get iteration history
// const vector<IterationData>& BisectionMethod::getHistory() const {
//     return history;
// }

// // Setters
// void BisectionMethod::setTolerance(double tol) {
//     tolerance = tol;
// }

// void BisectionMethod::setMaxIterations(int maxIter) {
//     maxIterations = maxIter;
// }

// void BisectionMethod::setInterval(double left, double right) {
//     a = left;
//     b = right;
// }

// // Getters
// double BisectionMethod::getTolerance() const {
//     return tolerance;
// }

// int BisectionMethod::getMaxIterations() const {
//     return maxIterations;
// }