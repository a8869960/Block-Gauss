//
// Created by varsem on 24.09.23.
//
#include "functions.h"
#include "gauss.h"

#include <cstring>

int gauss_func(int n,
               [[maybe_unused]]int m,
               double *A,
               [[maybe_unused]]double *B,
               double *x,
               double *Ahelp,
               double *inverseA,
               int *indi_m,
               int *indj_m,
               [[maybe_unused]]int *indi,
               [[maybe_unused]]int *indj,
               double* a)
{
    [[maybe_unused]]int k, l, bl, i, j;

    k = n / m; //how many blocks m*m
    l = n - k * m; //how long last block
    bl = (l != 0) ? k + 1 : k; //number of all blocks

    for(i = 0; i < bl; i++)
    {
        indi[i] = i;
        indj[i] = i;
    }

    memcpy(a, A, sizeof(double) * n * n);

    for(int step = 0; step < n; step++)
    {
//        if(matrixMax(a, step, n, m, k, l, indi, indj, indi_m, indj_m, block, inverse, help) == -1)
//            return -1;

        return 0;
    }

    if(inverseMatrix(A, inverseA, Ahelp, n, indi_m, indj_m) == -1)
        return -1;

    for(int i = 0; i < n; i++)
        x[i] = (i + 1) % 2;

    return 0;
}