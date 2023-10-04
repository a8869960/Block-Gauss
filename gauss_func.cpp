//
// Created by varsem on 24.09.23.
//
#include "functions.h"

int gauss_func(int n,
               [[maybe_unused]]int m,
               [[maybe_unused]]double *A,
               [[maybe_unused]]double *B,
               double *x,
               [[maybe_unused]]double *helper)
{
    for(int i = 0; i < n; i++)
        x[i] = (i + 1) % 2;

    if(helper) delete[] helper;

    return 0;
}

int matrixProduct(double* A1,
                  int n,
                  int m,
                  double *A2,
                  int r,
                  int s,
                  double *C)
{
    if(m != r)
    {
        cout << "Can't calculate the product of matrix." << endl;
        return -1;
    }

    for(int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= s; j++)
        {
            C[s * (i - 1) + j - 1] = 0;
            for(int k = 1; k <= m; k++)
            {
                C[s * (i - 1) + j - 1] += A1[m * (i - 1) + k - 1] * A2[s * (k - 1) + j - 1];
            }
        }

    }
    return 0;
}

int matrixSubtraction(double* A1,
                      int n,
                      int m,
                      double *A2,
                      int r,
                      int s,
                      double *C)
{
    if(n != r or m != s)
    {
        cout << "Can't do matrix subtraction." << endl;
        return -1;
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
            C[m * i + j] = A1[m * i + j] - A2[m * i + j];
    }

    return 0;
}