#ifndef EIGENVALUE_HPP
#define EIGENVALUE_HPP
#include "Matrix.hpp"
#include <vector>
#include <string>

class EigenValue : public Matrix {

public:
    // Constructs an n×n zero matrix ready to be filled.
    EigenValue(int n);
    // Copies from an existing Matrix object. Throws if not square.
    EigenValue(const Matrix& m);

    double powerMethod(std::vector<double>& eigenvec,
                       int maxIter = 10000,
                       double tol  = 1e-10) const;

    double inversePowerMethod(std::vector<double>& eigenvec,
                              int maxIter = 10000,
                              double tol  = 1e-10) const;

    std::vector<double> qrAlgorithm(int maxIter = 2000,
                                    double tol  = 1e-10) const;

    struct GershgorinDisk {
        int    row;
        double center;   // a_ii
        double radius;   // Σ_{j≠i} |a_ij|
        double lo() const { return center - radius; }
        double hi() const { return center + radius; }
    };
    std::vector<GershgorinDisk> gershgorinBounds() const;

    double conditionNumber() const;

    int getRows() const { return rows; }
    int getCols() const { return cols; }

private:
    // QR decomposition via modified Gram-Schmidt: A = Q·R
    void qrDecompose(const std::vector<std::vector<double>>& A,std::vector<std::vector<double>>& Q, std::vector<std::vector<double>>& R) const;

    // Gaussian elimination with partial pivoting: solves M·x = rhs
    std::vector<double> gaussSolve(std::vector<std::vector<double>> M,
                                   std::vector<double> rhs) const;

    // Matrix–vector multiply: returns A·v
    std::vector<double> matvec(const std::vector<std::vector<double>>& A,
                               const std::vector<double>& v) const;

    // Vector helpers
    double vecNorm(const std::vector<double>& v) const;
    void   vecNormalize(std::vector<double>& v) const;
    double vecDot(const std::vector<double>& a,
                  const std::vector<double>& b) const;

    // Copy internal data to a plain 2-D array (for algorithms)
    std::vector<std::vector<double>> toRaw() const;
};

#endif // EIGENVALUE_HPP