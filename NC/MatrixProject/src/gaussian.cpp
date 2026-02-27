#include "../include/gaussian.hpp"
#include <cmath>
#include <stdexcept>

GaussianElimination::GaussianElimination(bool usePivot)
{
    pivoting = usePivot;
}

vector<double> GaussianElimination::solve(Matrix &A)
{
    int n = A.getRows();

    //Forward Elimination
    for (int k = 0; k < n; k++)
    {
        // Partial Pivoting: swap in the row with the largest |value|
        if (pivoting)
        {
            int maxRow = k;
            double maxVal = fabs(A.get(k, k));

            for (int i = k + 1; i < n; i++)
            {
                if (fabs(A.get(i, k)) > maxVal)
                {
                    maxVal = fabs(A.get(i, k));
                    maxRow = i;
                }
            }
            A.swapRows(k, maxRow);
        }

        // Check for zero pivot (singular or near-singular matrix)
        if (fabs(A.get(k, k)) < 1e-12)
            throw runtime_error("Zero pivot encountered — system may have no unique solution.");

        // Eliminate all rows below row k
        for (int i = k + 1; i < n; i++)
        {
            double factor = A.get(i, k) / A.get(k, k);      

            for (int j = k; j < n + 1; j++)
                A.set(i, j, A.get(i, j) - factor * A.get(k, j)); 
        }
    }

    //Back Substitution
    vector<double> x(n);

    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = A.get(i, n); //compute RHS 

        for (int j = i + 1; j < n; j++)
            x[i] -= A.get(i, j) * x[j];  //substract known terms from RHS

        x[i] /= A.get(i, i);   //divide by pivot
    }  

    return x;
}
