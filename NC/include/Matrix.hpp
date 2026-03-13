#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>
#include <fstream>

class Matrix {

protected:
    int rows;
    int cols;
    std::vector<std::vector<double>> data;

public:
    Matrix();
    Matrix(int rows,int cols);
    Matrix(const Matrix &m); 

    void readFromFile(std::ifstream &fin);
    void displayToFile(std::ofstream &fout) const;

    Matrix operator+(const Matrix &m) const;
    Matrix operator-(const Matrix &m) const;
    Matrix operator*(const Matrix &m) const;

    double& operator()(int i,int j);//  for non-const objects
    double operator()(int i,int j) const;

    bool operator==(const Matrix &m) const;

    friend std::istream& operator>>(std::istream &in, Matrix &m);
    friend std::ostream& operator<<(std::ostream &out,const Matrix &m);

    bool isSquare() const;
    bool isSymmetric() const;
    bool isIdentity() const;
    bool isNull() const;
    bool isDiagonal() const;
    bool isDiagonallyDominant() const;
    bool isTranspose(const Matrix &m) const;

    void makeDiagonallyDominant();
    Matrix transpose() const;

    double determinant() const;
    Matrix inverse() const;
};

#endif
