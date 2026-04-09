#ifndef SLE_HPP
#define SLE_HPP
#include "Matrix.hpp"
#include <vector>

class SLE : public Matrix {
public:
    SLE(int r,int c);
    virtual ~SLE();
     virtual std::vector<double> solve() = 0;

 protected:
    void prepareForIteration(const std::string& methodName);
};

#endif