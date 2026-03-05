#include "../include/LU.hpp"
#include <cmath>
#include <vector>

LU::LU(int r, int c) : SLE(r, c) {}

void LU::solve(std::ofstream& outFile) {
    int n = rows;
   
    // STEP 1: Doolittle LU Decomposition
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {

        // Fill row i of U (upper triangular) 
        for (int k = i; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j)
                sum += L[i][j] * U[j][k];
            U[i][k] = data[i][k] - sum;
        }

        // A zero here means the matrix is singular (no unique solution).
        if (std::abs(U[i][i]) < 1e-12) {
            outFile << "Zero (or near-zero) pivot U[" << i << "][" << i
                    << "]. Matrix is singular. Cannot solve.\n";
            return;
        }

        //Fill column i of L (lower triangular)
        for (int k = i; k < n; ++k) {
            if (i == k) {
                L[i][i] = 1.0; // Diagonal of L is always 1 in Doolittle's method
            } else {
                double sum = 0.0;
                for (int j = 0; j < i; ++j)
                    sum += L[k][j] * U[j][i];
                L[k][i] = (data[k][i] - sum) / U[i][i];
            }
        }
    }

    // STEP 2: Forward Substitution — solve Ly = b
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = data[i][n]; // b_i (last column of augmented matrix)
        for (int j = 0; j < i; ++j)
            y[i] -= L[i][j] * y[j];
        // No division: L[i][i] == 1.0 always
    }

    // STEP 3: Backward Substitution — solve Ux = y
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= U[i][j] * x[j];
        x[i] /= U[i][i]; // Safe: zero pivot already checked above
    }

    outFile << "Solution using LU Decomposition (Doolittle):\n";
    for (int i = 0; i < n; ++i)
        outFile << "x" << i + 1 << " = " << x[i] << "\n";
}