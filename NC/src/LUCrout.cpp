#include "../include/LUCrout.hpp"

LUCrout::LUCrout(int r,int c):SLE(r,c){}

std::vector<double> LUCrout::solve()
{
    int n=rows;

    std::vector<std::vector<double>> L(n,std::vector<double>(n,0));
    std::vector<std::vector<double>> U(n,std::vector<double>(n,0));

    for(int j=0;j<n;j++)
    {
        U[j][j]=1;

        for(int i=j;i<n;i++)
        {
            double sum=0;

            for(int k=0;k<j;k++)
                sum+=L[i][k]*U[k][j];

            L[i][j]=data[i][j]-sum;
        }

        for(int i=j+1;i<n;i++)
        {
            double sum=0;

            for(int k=0;k<j;k++)
                sum+=L[j][k]*U[k][i];

            U[j][i]=(data[j][i]-sum)/L[j][j];
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
            x[i]-=U[i][j]*x[j];
    }

    return x;
}