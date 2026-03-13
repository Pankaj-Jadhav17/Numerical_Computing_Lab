#include "../include/GaussElimination.hpp"

GaussElimination::GaussElimination(int r,int c)
: SLE(r,c)
{
}

std::vector<double> GaussElimination::solve()
{
    int n = rows;
    for(int k=0;k<n-1;k++)
    {
        for(int i=k+1;i<n;i++)
        {
            double factor = data[i][k]/data[k][k];
            for(int j=k;j<cols;j++)
            {
                data[i][j] -= factor * data[k][j];
            }
        }
    }

    std::vector<double> x(n);

    for(int i=n-1;i>=0;i--)
    {
        x[i] = data[i][n];
        for(int j=i+1;j<n;j++)
        {
            x[i] -= data[i][j] * x[j];
        }
        x[i] /= data[i][i];
    }
    return x;
}