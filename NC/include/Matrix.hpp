#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <fstream>

class Matrix {
protected:
    int rows;
    int cols;
    std::vector<std::vector<double>> data; 

public:
    Matrix(int r, int c); 
    virtual ~Matrix() = default; 
    // Read matrix values row by row from a file
    void readFromFile(std::ifstream& inFile);

    // Write matrix values to a file (useful for debugging)
    void writeMatrix(std::ofstream& outFile) const;

    int getRows() const;
    int getCols() const;
};

#endif 