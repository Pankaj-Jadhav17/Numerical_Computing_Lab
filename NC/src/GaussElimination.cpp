#include "../include/GaussElimination.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

GaussElimination::GaussElimination(int r, int c, bool pivot)
    : SLE(r, c), usePivoting(pivot) {}

void GaussElimination::solve(std::ofstream& outFile) {
    int n = rows;
    // data is the augmented matrix [A | b]:
    // STEP 1: Forward Elimination
    for (int i = 0; i < n - 1; ++i) {

        if (usePivoting) {
            int maxRow = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(data[k][i]) > std::abs(data[maxRow][i]))
                    maxRow = k;
            }
            if (maxRow != i)
                std::swap(data[i], data[maxRow]);
        }

        // Check for zero pivot — means the system has no unique solution
        if (std::abs(data[i][i]) < 1e-12) {
            outFile << "Zero (or near-zero) pivot at row " << i
                    << ". Matrix may be singular. Cannot solve.\n";
            return;
        }

        // Eliminate all entries below the pivot in column i
        for (int j = i + 1; j < n; ++j) {
            double ratio = data[j][i] / data[i][i];
            for (int k = i; k < cols; ++k)  // start k at i (no need to update already-zero entries)
                data[j][k] -= ratio * data[i][k];
        }
    }

    // Final pivot check for the last row
    if (std::abs(data[n-1][n-1]) < 1e-12) {
        outFile << "Zero pivot at the last row. Matrix is singular. Cannot solve.\n";
        return;
    }

    // STEP 2: Back Substitution
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = data[i][cols - 1]; // right-hand side value (b_i)
        for (int j = i + 1; j < n; ++j)
            x[i] -= data[i][j] * x[j]; // subtract already-known x values
        x[i] /= data[i][i];
    }

    outFile << "Solution using Gauss Elimination"
            << (usePivoting ? " (with partial pivoting)" : " (no pivoting)") << ":\n";
    for (int i = 0; i < n; ++i)
        outFile << "x" << i + 1 << " = " << x[i] << "\n";
}