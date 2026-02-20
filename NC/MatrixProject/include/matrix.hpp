#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;

// Abstract Parent Class
class MatrixBase
{
protected:
    int rows, cols;

public:
    MatrixBase(); 
    MatrixBase(int r, int c);

    //Virtual destructor
    virtual ~MatrixBase();

    //Pure virtual functions → abstract class
    virtual void read() = 0;
    virtual void display() const = 0;

    virtual MatrixBase* add(const MatrixBase &m) const = 0;
    virtual MatrixBase* subtract(const MatrixBase &m) const = 0;

    virtual void gaussianElimination() = 0;
    virtual void backSubstitution() = 0;
};

//Child Class
class Matrix : public MatrixBase
{
private:
    vector<vector<double>> mat; // 2D vector to store matrix elements

public:
    // Constructors
    Matrix();
    Matrix(int r, int c);

    // Copy constructor
    Matrix(const Matrix &other);

    // Destructor
    ~Matrix();

    // Overriding
    void read() override;
    void display() const override;

    MatrixBase* add(const MatrixBase &m) const override;
    MatrixBase* subtract(const MatrixBase &m) const override;

    void gaussianElimination() override;
    void backSubstitution() override;
};

#endif






































// #ifndef MATRIX_HPP
// #define MATRIX_HPP

// #include <iostream>
// #include <vector>
// #include <stdexcept>
// using namespace std;

// // 🔹 Parent class
// class MatrixBase
// {
// protected:
//     int rows, cols;

// public:
//     MatrixBase() : rows(0), cols(0) {}

//     // 🔹 Virtual constructor concept (factory)
//     virtual MatrixBase* create(int r, int c) = 0;

//     // Virtual functions
//     virtual void read() = 0;
//     virtual void display() = 0;

//     virtual MatrixBase* add(MatrixBase &m) = 0;
//     virtual MatrixBase* subtract(MatrixBase &m) = 0;

//     virtual void gaussianElimination() = 0;
//     virtual void backSubstitution() = 0;

//     // 🔹 Virtual destructor
//     virtual ~MatrixBase()
//     {
//         cout << "Base class destructor\n";
//     }
// };

// // 🔹 Child class
// class Matrix : public MatrixBase
// {
// private:
//     vector<vector<double>> mat;

// public:
//     Matrix();
//     Matrix(int r, int c);

//     // overriding factory
//     MatrixBase* create(int r, int c) override;

//     void read() override;
//     void display() override;

//     MatrixBase* add(MatrixBase &m) override;
//     MatrixBase* subtract(MatrixBase &m) override;

//     void gaussianElimination() override;
//     void backSubstitution() override;

//     // destructor
//     ~Matrix()
//     {
//         cout << "Matrix destructor\n";
//     }
// };

// #endif
