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
    double *Ahelp = new double[n*n], *inverseA = new double[n*n];
    int *indi = new int[n], *indj = new int[n];

//    cout << "--------invA" << endl;
    inverseMatrix(A, Ahelp, inverseA, n, indi, indj);
//    matrixOutput(inverseA, n, n, n);
//    cout << "--------Ahelp" << endl;
//    matrixOutput(Ahelp, n, n, n);
//    cout << "--------A" << endl;
//    matrixOutput(A, n, n, n);

    for(int i = 0; i < n; i++)
        x[i] = (i + 1) % 2;

    if(helper) delete[] helper;

    return 0;
}