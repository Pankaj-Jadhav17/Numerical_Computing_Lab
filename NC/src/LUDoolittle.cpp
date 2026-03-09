#include "../include/LUDoolittle.hpp"

LUDoolittle::LUDoolittle(int r,int c):SLE(r,c){}

std::vector<double> LUDoolittle::solve()
{
    int n=rows;

    std::vector<std::vector<double>> L(n,std::vector<double>(n,0));
    std::vector<std::vector<double>> U(n,std::vector<double>(n,0));

    for(int i=0;i<n;i++)
    {
        for(int k=i;k<n;k++)
        {
            double sum=0;

            for(int j=0;j<i;j++)
                sum+=L[i][j]*U[j][k];

            U[i][k]=data[i][k]-sum;
        }

        for(int k=i;k<n;k++)
        {
            if(i==k)
                L[i][i]=1;
            else
            {
                double sum=0;

                for(int j=0;j<i;j++)
                    sum+=L[k][j]*U[j][i];

                L[k][i]=(data[k][i]-sum)/U[i][i];
            }
        }
    }

    std::vector<double> y(n);

    for(int i=0;i<n;i++)
    {
        y[i]=data[i][n];

        for(int j=0;j<i;j++)
            y[i]-=L[i][j]*y[j];
    }

    std::vector<double> x(n);

    for(int i=n-1;i>=0;i--)
    {
        x[i]=y[i];

        for(int j=i+1;j<n;j++)
            x[i]-=U[i][j]*x[j];

        x[i]/=U[i][i];
    }

    return x;
}