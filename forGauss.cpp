//
// Created by varsem on 09.10.23.
//
#include "gauss.h"

#include <cstring>

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
              double *block_inv,
              double *block_h)
{
    double min = 1.7976931348623158e+308, norm;
    int imax = step, jmax = step, count = 0;

    for(int i = step; i < k; i++)
        for(int j = step; j < k; j++) {
            get_block(A, block, indi[i], indj[j], n, m, k, l);

            if (inverseMatrix(block, block_inv, block_h, m, indi_m, indj_m) == 0)
            {
                norm = matrixNorm(block_inv, m);
//                cout << "NORM " << norm << endl;

                if (norm < min)
                {
                    min = norm;

                    imax = i;
                    jmax = j;
//                    cout << "MAX " << imax << jmax << endl;
                }
            } else
                count++;
        }

    if(count == (k - step) * (k - step))
        return -1;

    int helper;

    helper = indi[imax];
    indi[imax] = indi[step];
    indi[step] = helper;

    helper = indj[jmax];
    indj[jmax] = indj[step];
    indj[step] = helper;
    
//    cout << "MAX " << imax << jmax << endl;

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
        int i,
        int j,
        int n,
        int m,
        int k,
        int l)
{
    int block_m = (i == k ? l : m), block_l = (j == k ? l : m);

//    cout << " i " << i << " " << "j " << j << " " << " block_ m = " << block_m << " " << "block_l = " << block_l << endl;
    int r, s;
    int a = i * n * m + j * m; //number of first element of the block

    for(r = 0; r < block_m; r++)
    {
        for(s = 0; s < block_l; s++)
            block[r * m + s] = A[a + r * n + s];
        for(s = block_l; s < m; s++)
            block[r * m + s] = 0;
    }


    for(r = block_m; r < m; r++)
        for(s = 0; s < m; s++)
            block[r * m + s] = 0;
}

void put_block(
        double* A,
        double* block,
        int i,
        int j,
        int n,
        int m,
        int k,
        int l)
{
    int block_m = (i == k ? l : m), block_l = (j == k ? l : m);

    int r, s;
    int a = i * n * m + j * m; //number of first element of the block

    for(r = 0; r < block_m; r++)
    {
        for(s = 0; s < block_l; s++)
        {
            A[a + r * n + s] = block[r * block_l + s];
        }
    }
}

void E(double* block, int m)
{
    memset(block, 0, sizeof(double) * m * m);

    for(int i = 0; i < m; i++)
        block[i * m + i] = 1;
}

void get_block_b( double *B, double *block, int i, int m, int k, int l)
{
    int block_m = (i == k ? l : m);

    int r;
    int b = i * m; //number of first element of the block

    for(r = 0; r < block_m; r++)
    {
        block[r] = B[b + r];
    }

    for(r = block_m; r < m; r++)
        block[r] = 0;
}

void put_block_b( double *B, double *block, int i, int m, int k, int l)
{
    int block_m = (i == k ? l : m);

    int r;
    int b = i * m; //number of first element of the block

    for(r = 0; r < block_m; r++)
    {
         B[b + r] = block[r];
    }
}