#include "../include/GaussSeidel.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

GaussSeidel::GaussSeidel(int r, int c) : SLE(r, c) {}

std::vector<double> GaussSeidel::solve()
{
    int n = rows;
        prepareForIteration("Gauss-Seidel");

    const int    maxIter   = 10000;  // enough headroom for large/slow systems
    const double tolerance = 1e-10;  // tight tolerance for accuracy

    std::vector<double> x(n, 0.0);      // initial guess: all zeros
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