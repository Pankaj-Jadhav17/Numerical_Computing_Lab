#include "../include/LUCholesky.hpp"
#include <cmath>

LUCholesky::LUCholesky(int r,int c):SLE(r,c){}

std::vector<double> LUCholesky::solve()
{
    int n=rows;

    std::vector<std::vector<double>> L(n,std::vector<double>(n,0));

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<=i;j++)
        {
            double sum=0;

            for(int k=0;k<j;k++)
                sum+=L[i][k]*L[j][k];

            if(i==j)
                L[i][j]=sqrt(data[i][i]-sum);
            else
                L[i][j]=(data[i][j]-sum)/L[j][j];
        }
    }

    std::vector<double> y(n);

    for(int i=0;i<n;i++)
    {
        y[i]=data[i][n];

        for(int j=0;j<i;j++)
            y[i]-=L[i][j]*y[j];

        y[i]/=L[i][i];
    }

    std::vector<double> x(n);

    for(int i=n-1;i>=0;i--)
    {
        x[i]=y[i];

        for(int j=i+1;j<n;j++)
            x[i]-=L[j][i]*x[j];

        x[i]/=L[i][i];
    }

    return x;
}