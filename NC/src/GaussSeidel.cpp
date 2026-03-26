#include "../include/GaussSeidel.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

GaussSeidel::GaussSeidel(int r, int c) : SLE(r, c) {}

std::vector<double> GaussSeidel::solve()
{
    int n = rows;

    // ── Step 1: Partial Pivoting — bring largest |value| to diagonal ──────
    // Improves diagonal dominance and convergence for large matrices
    for (int col = 0; col < n; col++)
    {
        int maxRow = col;
        for (int row = col + 1; row < n; row++)
            if (std::fabs(data[row][col]) > std::fabs(data[maxRow][col]))
                maxRow = row;

        if (maxRow != col)
            std::swap(data[col], data[maxRow]);
    }

    // ── Step 2: Check diagonal dominance (warn if not met) ───────────────
    bool dominant = true;
    for (int i = 0; i < n; i++)
    {
        double diag = std::fabs(data[i][i]);
        double sum  = 0.0;
        for (int j = 0; j < n; j++)
            if (j != i) sum += std::fabs(data[i][j]);

        if (diag < sum) { dominant = false; break; }
    }
    if (!dominant)
        std::cout << "[Warning] Matrix is not strictly diagonally dominant."
                  << " Gauss-Seidel may not converge!\n";

    std::vector<double> x(n, 0.0);   // initial guess: all zeros

    const int    maxIter   = 10000;  // enough headroom for large/slow systems
    const double tolerance = 1e-10;  // tight tolerance for accuracy

    for (int iter = 0; iter < maxIter; iter++)
    {
        double maxChange = 0.0;

        for (int i = 0; i < n; i++)
        {
            double sigma = data[i][n];   // start with b[i]
            for (int j = 0; j < n; j++)
                if (j != i)
                    sigma -= data[i][j] * x[j];

            double x_new = sigma / data[i][i];

            // Track the largest change for convergence check
            maxChange = std::max(maxChange, std::fabs(x_new - x[i]));
            x[i] = x_new;   // update immediately (in-place)
        }

        if (maxChange < tolerance)
        {
            std::cout << "(Converged after " << iter + 1 << " iterations)\n";
            return x;
        }
    }
    std::cout << "[Warning] Gauss-Seidel did not converge in "
              << maxIter << " iterations. Returning best estimate.\n";
    return x;
}