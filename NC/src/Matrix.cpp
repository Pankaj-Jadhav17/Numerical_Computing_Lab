#include "../include/Matrix.hpp"

Matrix::Matrix(int r, int c) : rows(r), cols(c) {
    // Initialize a 2D grid of zeros with dimensions rows x cols
    data.resize(rows, std::vector<double>(cols, 0.0));
}

void Matrix::readFromFile(std::ifstream& inFile) {
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            inFile >> data[i][j];
}

void Matrix::writeMatrix(std::ofstream& outFile) const { 
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j)
            outFile << data[i][j] << " ";
        outFile << "\n";
    }
}

int Matrix::getRows() const { return rows; }
int Matrix::getCols() const { return cols; }