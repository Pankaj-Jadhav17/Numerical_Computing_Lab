#include "../include/SLE.hpp"
#include <cmath>    

SLE::SLE(int r, int c) : Matrix(r, c) {}
SLE::~SLE() {}

void SLE::prepareForIteration(const std::string& methodName)
{
    makeDiagonallyDominant();   // partial pivoting
 
    if (!isDiagonallyDominant())  // dominance check 
        std::cout<< "[Warning] Matrix is not diagonally dominant."<< methodName <<"may not converge!\n";
}
 