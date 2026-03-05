#ifndef SLE_H
#define SLE_H

#include "Matrix.hpp"

// Any solver (Gauss, LU, etc.) must inherit this class and implement solve().
class SLE : public Matrix {
public:
    SLE(int r, int c) : Matrix(r, c) {}  
        virtual ~SLE() = default; 

    // Pure virtual: each subclass must provide its own solve() method
    virtual void solve(std::ofstream& outFile) = 0;
};

#endif 