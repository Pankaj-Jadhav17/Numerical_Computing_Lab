#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>

using namespace std;

class Matrix
{
protected:
    int rows, cols;
    vector<vector<double>> data;

public:
    Matrix(int r = 0, int c = 0);

    int getRows() const;
    int getCols() const;

    double get(int i, int j) const;
    void set(int i, int j, double value);

    void swapRows(int r1, int r2);

    // File utilities
    void readFromFile(const string &file);
    void writeSolution(const string &file, const vector<double> &x);
};

#endif

















// #ifndef MATRIX_HPP
// #define MATRIX_HPP

// #include <vector>
// #include <stdexcept>
// using namespace std;

// class Matrix
// {
// private: 
//     int rows, cols; 
//     vector<vector<double>> data;

// public:
//     Matrix(int r = 0, int c = 0);

//     int getRows() const;
//     int getCols() const;

//     double get(int i, int j) const;
//     void set(int i, int j, double value);
//     void swapRows(int r1, int r2);
// };

// #endif







// #ifndef MATRIX_HPP
// #define MATRIX_HPP

// #include <iostream>
// #include <vector>
// #include <stdexcept>
// using namespace std;

// // Abstract Parent Class
// class MatrixBase
// {
// protected:
//     int rows, cols;

// public:
//     MatrixBase(); 
//     MatrixBase(int r, int c);

//     //Virtual destructor
//     virtual ~MatrixBase();

//     //Pure virtual functions → abstract class
//     virtual void read() = 0;
//     virtual void display() const = 0;

//     virtual MatrixBase* add(const MatrixBase &m) const = 0;
//     virtual MatrixBase* subtract(const MatrixBase &m) const = 0;

//     virtual void gaussianElimination() = 0;
//     virtual void backSubstitution() = 0;
// };

// //Child Class
// class Matrix : public MatrixBase
// {
// private:
//     vector<vector<double>> mat; // 2D vector to store matrix elements

// public:
//     // Constructors
//     Matrix();
//     Matrix(int r, int c);

//     // Copy constructor
//     Matrix(const Matrix &other);

//     // Destructor
//     ~Matrix();

//     // Overriding
//     void read() override;
//     void display() const override;

//     MatrixBase* add(const MatrixBase &m) const override;
//     MatrixBase* subtract(const MatrixBase &m) const override;

//     void gaussianElimination() override;
//     void backSubstitution() override;
// };

// #endif






























