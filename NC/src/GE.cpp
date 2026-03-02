#include "../include/GE.hpp"
#include <cmath>
#include <stdexcept>

GE::GE(int n, bool usePivot) : SLE(n)
{
    pivoting = usePivot;
}

vector<double> GE::solve()
{
    int n = rows;

    for (int k = 0; k < n; k++)
    {
        if (pivoting)
        {
            int maxRow = k;
            double maxVal = fabs(get(k, k));

            for (int i = k + 1; i < n; i++)
            {
                if (fabs(get(i, k)) > maxVal)
                {
                    maxVal = fabs(get(i, k));
                    maxRow = i;
                }
            }
            swapRows(k, maxRow);
        }

        if (fabs(get(k, k)) < 1e-12)
            throw runtime_error("Zero pivot");

        for (int i = k + 1; i < n; i++)
        {
            double factor = get(i, k) / get(k, k);

            for (int j = k; j < n + 1; j++)
                set(i, j, get(i, j) - factor * get(k, j));
        }
    }

    vector<double> x(n);

    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = get(i, n);

        for (int j = i + 1; j < n; j++)
            x[i] -= get(i, j) * x[j];

        x[i] /= get(i, i);
    }

    return x;
}