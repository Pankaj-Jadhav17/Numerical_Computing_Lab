#include "../include/GaussJacobi.hpp"
#include <cmath>
#include <algorithm>  // for swap

GaussJacobi::GaussJacobi(int r, int c) : SLE(r, c) {}

std::vector<double> GaussJacobi::solve()
{
    int n = rows;

    //Step 1: Reorder rows so the largest value is on the diagonal
    // This helps Gauss-Jacobi converge (called partial pivoting)
    for (int col = 0; col < n; col++)
    {
        int maxRow = col;
        for (int row = col + 1; row < n; row++)
            if (std::fabs(data[row][col]) > std::fabs(data[maxRow][col]))
                maxRow = row;

        std::swap(data[col], data[maxRow]);   // swap entire rows
    }

 // Step 2: Check diagonal dominance (warn user if not met)
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
        std::cout << "[Warning] Matrix is not diagonally dominant."
                  << " Jacobi may not converge!\n";

    // ── Step 3: Iterate until convergence ────────────────────────────
    std::vector<double> x(n, 0.0);      // start with x = [0, 0, ...]
    std::vector<double> x_new(n, 0.0);

    int    maxIter   = 10000;
    double tolerance = 1e-10;

    for (int iter = 0; iter < maxIter; iter++)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = data[i][n];           // start with b[i]
            for (int j = 0; j < n; j++)
                if (j != i)
                    sum -= data[i][j] * x[j]; // subtract using OLD x values

            x_new[i] = sum / data[i][i];
        }

        // Check if values stopped changing (convergence)
        double maxChange = 0.0;
        for (int i = 0; i < n; i++)
            maxChange = std::max(maxChange, std::fabs(x_new[i] - x[i]));

        x = x_new;

        if (maxChange < tolerance)
        {
            std::cout << "(Converged after " << iter + 1 << " iterations)\n";
            return x;
        }
    }

    std::cout << "[Warning] Gauss-Jacobi did not converge in "
              << maxIter << " iterations. Returning best estimate.\n";
    return x;
}