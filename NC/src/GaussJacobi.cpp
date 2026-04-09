#include "../include/GaussJacobi.hpp"
#include <cmath>
#include <algorithm>  // for swap

GaussJacobi::GaussJacobi(int r, int c) : SLE(r, c) {}

std::vector<double> GaussJacobi::solve()
{
    int n = rows;
    prepareForIteration("Jacobi");
    std::vector<double> x(n, 0.0);      // initial guess: all zeros
    std::vector<double> x_new(n, 0.0);  

    const int    maxIter   = 10000;  // enough headroom for large/slow systems
    const double tolerance = 1e-10;  // tight tolerance for accuracy


    for (int iter = 0; iter < maxIter; iter++)
    {
        for(int i = 0; i < n; i++)
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