#ifndef BM_HPP
#define BM_HPP
#include <vector>

struct IterationData {
    int iteration;
    double a;
    double b;
    double midpoint;
    double f_midpoint;
    double error;
    double interval_length;
};
class BisectionMethod {
private:
    double a, b;
    double tolerance;
    int k;
    std::vector<IterationData> history;
    double f(double x);
public:
    BisectionMethod(double left, double right, int tolerance_power);
    bool validateInterval();
    double findRoot();
    void displayIterationTable();
    double getTolerance() const;
    int getActualIterations() const;
};

#endif














/*
================================================================================
    BM.hpp - Header File for Bisection Method
    Function: f(x) = 4x³ - 3x
================================================================================
*/

// #ifndef BM_HPP
// #define BM_HPP

// #include <vector>
// #include <string>

// // Structure to store iteration data
// struct IterationData {
//     int iteration;              // Iteration number
//     double a;                   // Left endpoint
//     double b;                   // Right endpoint
//     double midpoint;            // Midpoint value
//     double f_midpoint;          // Function value at midpoint
//     double error;               // |m_current - m_previous|
//     double interval_length;     // b - a
// };

// // BisectionMethod class for f(x) = 4x³ - 3x
// class BisectionMethod {
// private:
//     double a;                                // Left endpoint
//     double b;                                // Right endpoint
//     double tolerance;                        // Stopping criterion
//     int k;                                   // Tolerance power (tolerance = 10^-k)
//     std::vector<IterationData> history;     // Store iteration history
    
//     // The function f(x) = 4x³ - 3x
//     double f(double x);
    
// public:
//     // Constructor
//     BisectionMethod(double left, double right, int tolerance_power);
    
//     // Step 1: Find the root using bisection method
//     double findRoot();
    
//     // Step 2: Calculate theoretical number of iterations
//     int calculateTheoreticalIterations();
    
//     // Step 3: Validate that interval brackets a root
//     bool validateInterval();
    
//     // Step 4: Get iteration history (shows interval at each iteration)
//     const std::vector<IterationData>& getIterationHistory() const;
    
//     // Step 5: Display iteration table showing interval reduction
//     void displayIterationTable();
    
//     // Step 6: Display final results including tolerance
//     void displayResults(double root, int actualIterations, int theoreticalIterations);
    
//     // Get tolerance value
//     double getTolerance() const;
    
//     // Get number of actual iterations
//     int getActualIterations() const;
    
//     // Evaluate the function f(x) = 4x³ - 3x
//     double evaluateFunction(double x) const;
// };

// #endif // BM_HPP





























// #ifndef BM_HPP
// #define BM_HPP

// #include <functional>
// #include <vector>
// #include <string>

// // Structure to store iteration data
// struct IterationData {
//     int iteration;
//     double midpoint;
//     double error;
//     double intervalLeft;
//     double intervalRight;
// };

// // BisectionMethod class for finding roots of equations
// class BisectionMethod {
// private:
//     std::function<double(double)> function;  // The equation to solve
//     double a;                                 // Left endpoint
//     double b;                                 // Right endpoint
//     double tolerance;                         // Stopping criterion
//     int maxIterations;                        // Maximum iterations allowed
//     std::vector<IterationData> history;      // Store iteration history
    
// public:
//     // Constructor
//     BisectionMethod(std::function<double(double)> func, 
//                     double left, 
//                     double right, 
//                     double tol = 1e-6,
//                     int maxIter = 1000);
    
//     // Main method to find the root
//     double findRoot();
    
//     // Calculate theoretical minimum iterations
//     int getTheoreticalIterations() const;
    
//     // Validate initial interval
//     bool validateInterval() const;
    
//     // Get iteration history
//     const std::vector<IterationData>& getHistory() const;
    
//     // Display results
//     void displayResults(double root, int actualIterations) const;
    
//     // Display iteration table
//     void displayIterationTable() const;
    
//     // Setters for parameters
//     void setTolerance(double tol);
//     void setMaxIterations(int maxIter);
//     void setInterval(double left, double right);
    
//     // Getters
//     double getTolerance() const;
//     int getMaxIterations() const;
// };

// #endif // BM_HPP