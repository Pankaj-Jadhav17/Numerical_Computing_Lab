#include "../include/GaussElimination.hpp"
#include <cmath>
#include <stdexcept>

GaussElimination::GaussElimination(int r,int c)
: SLE(r,c)
{
}

std::vector<double> GaussElimination::solve()
{
    int n = rows;

    // ── Forward Elimination with Partial Pivoting ──────────────
    for(int k = 0; k < n-1; k++)
    {
        // Step 1: find the row with the largest |pivot| in column k
        int maxRow = k;
        for(int i = k+1; i < n; i++)
            if(std::fabs(data[i][k]) > std::fabs(data[maxRow][k]))
                maxRow = i;

        // Step 2: swap rows k and maxRow (entire augmented row)
        if(maxRow != k)
            std::swap(data[k], data[maxRow]);

        // Step 3: guard against a zero (or near-zero) pivot
        if(std::fabs(data[k][k]) < 1e-12)
            throw std::runtime_error(
                "GaussElimination: singular or near-singular matrix at pivot " +
                std::to_string(k));

        // Step 4: eliminate entries below the pivot
        for(int i = k+1; i < n; i++)
        {
            double factor = data[i][k] / data[k][k];
            for(int j = k; j < cols; j++)
                data[i][j] -= factor * data[k][j];
        }
    }

    // ── Back Substitution ──────────────────────────────────────
    std::vector<double> x(n);
    for(int i = n-1; i >= 0; i--)
    {
        x[i] = data[i][n];          // last column = b (RHS)
        for(int j = i+1; j < n; j++)
            x[i] -= data[i][j] * x[j];
        x[i] /= data[i][i];
    }
    return x;
}