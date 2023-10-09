//
// Created by varsem on 09.10.23.
//
#include "gauss.h"

int matrixMax(double *A,
              int step,
              int n,
              int m,
              int k,
              int l,
              int *indi,
              int *indj,
              int *indi_m,
              int *indj_m,
              double* block,
              double *inverse,
              double *help)
{
    double min = 1.7976931348623158e+308, norm;
    int imax = step, jmax = step, count = 0;

    for(int i = step; i < n; i++)
        for(int j = step; j < n; j++) {
            get_block(A, block, indi[i], indj[j], n, m, k, l);

            if (inverseMatrix(block, inverse, help, m, indi_m, indj_m) == 0)
            {
                norm = matrixNorm(inverse, m);

                if (norm < min)
                {
                    min = norm;

                    imax = i;
                    jmax = j;
                }
            } else
                count++;
        }

    if(count == step * step)
        return -1;

    int helper;

    helper = indi[imax];
    indi[imax] = indi[step];
    indi[step] = helper;

    helper = indj[jmax];
    indj[jmax] = indj[step];
    indj[step] = helper;

    return 0;
}

double matrixNorm(double *A, int n)
{
    double norm = 0, helper = 0;

    for(int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
            helper += abs(A[i * n + j]);

        if(helper > norm)
            norm = helper;

        helper = 0;
    }

    return norm;
}

void get_block(
        double* A,
        double* block,
        size_t i,
        size_t j,
        size_t n,
        size_t m,
        size_t k,
        size_t l)
{
    int block_m = (i > k ? l : m), block_l = (j > k ? l : m);

    int r, s;
    int a = i * n * m + j * m; //number of first element of the block

    for(r = 0; r < block_m; r++)
    {
        for(s = 0; s < block_l; s++)
        {
            block[r * block_m + s] = A[a + r * n + s];
        }
    }
}

void put_block(
        double* A,
        double* block,
        size_t i,
        size_t j,
        size_t n,
        size_t m,
        size_t k,
        size_t l)
{
    int block_m = (i > k ? l : m), block_l = (j > k ? l : m);

    int r, s;
    int a = i * n * m + j * m; //number of first element of the block

    for(r = 0; r < block_m; r++)
    {
        for(s = 0; s < block_l; s++)
        {
            A[a + r * n + s] = block[r * block_m + s];
        }
    }
}