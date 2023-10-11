//
// Created by varsem on 24.09.23.
//
#include "functions.h"
#include "gauss.h"

#include <cstring>

int gauss_func(int n,
               int m,
               double *A,
               double *B,
               double *x,
               [[maybe_unused]]double *Ahelp,
               [[maybe_unused]]double *inverseA,
               int *indi_m,
               int *indj_m,
               int *indi,
               int *indj,
               double* a,
               double *b,
               double *block,
               double *block_inv,
               double *block_h,
               [[maybe_unused]]double *block_ii)
{
    [[maybe_unused]]int k, l, bl, i, j, step;

    k = n / m; //how many blocks m*m
    l = n - k * m; //how long last block
    bl = (l != 0) ? k + 1 : k; //number of all blocks

    for(i = 0; i < bl; i++)
    {
        indi[i] = i;
        indj[i] = i;
    }

    memcpy(a, A, sizeof(double) * n * n);
    memcpy(b, B, sizeof(double) * n);

    //Прямой ход
    for(step = 0; step < k; step++)
    {
        if(matrixMax(a, step, n, m, k, l, indi, indj, indi_m, indj_m, block, block_inv, block_h) == -1)
        {
            cout << "Can't find block max." << endl;
            return -1;
        }

        get_block(a, block, indi[step], indj[step], n, m, k, l);
        inverseMatrix(block, block_inv, block_h, m, indi_m, indj_m);

        cout << "BLOCK II" << endl;
        matrixOutput(block_inv, m, m, 3);

        //Деление первой строчки
        E(block, m);
        put_block(a, block, indi[step], indj[step], n, m, k, l);
        for(j = step + 1; j < bl; j++)
        {
            get_block(a, block, indi[step], indj[j], n, m, k, l);
//            cout << "BLOCK" << endl;
//            matrixOutput(block, m, m, 3);
            unit(block_inv, block, block_h, m, m);
//            cout << "BLOCK H" << endl;
//            matrixOutput(block_h, m, m, 3);
            put_block(a, block_h, indi[step], indj[j], n, m, k, l);
        }
        get_block_b(b, block, indi[step], m, k, l);
        matrixProduct(block_inv, m, m, block, m, 1, block_h);
        put_block_b(b, block_h, indi[step], m, k, l);

        //Зануление столбца
        for(i = step + 1; i < bl; i++)
        {
            get_block(a, block, indi[i], indj[step], n, m, k ,l);
            for(j = step + 1; j < bl; j++)
            {
                get_block(a, block_h, indi[step], indj[j], n, m, k, l);
                unit(block, block_h, block_inv, n, m);
                get_block(a, block_h, indi[i], indj[j], n, m, k, l);
                matrixSubtraction(block_h, m, m, block_inv, m, m, block_h);
                put_block(a, block_h, indi[i], indj[j], n, m, k, l);
            }
            get_block_b(b, block_h, indi[step], m, k, l);
            matrixProduct(block, m, m, block_h, m, 1, block_inv);
            get_block_b(b, block_h, indi[i], m, k, l);
            matrixSubtraction(block_h, 1, m, block_inv, 1, m, block_h);
            put_block_b(b, block_h, indi[i], m, k, l);
        }
        memset(block, 0, sizeof(double) * m * m);
        for(i = step + 1; i < bl; i++)
            put_block(a, block, indi[i], indj[step], n, m, k, l);

        cout << "PR HOD " << endl;
        matrixOutput(a, n, n, 3);
        cout << endl;
    }

    //Обратный ход
    for(step = bl - 1; step >= 0; step--)
    {
        get_block(a, block, indi[step], indj[step], n, m, k, l);
        inverseMatrix(block, block_inv, block_h, m, indi_m, indj_m);
        E(block, m);
        put_block(a, block, indi[step], indj[step], n, m, k, l);

        get_block_b(b, block, indi[step], m, k, l);
        matrixProduct(block_inv, m, m, block, m, 1, block_h);
        put_block_b(b, block_h, indi[step], m, k, l);

        for(i = bl - 2; i >= 0; i--)
        {
            get_block(a, block_h, indi[i], indj[step], n, m, k, l);
            matrixProduct(block_h, m, m, block, m, 1, block_inv);
            get_block_b(b, block, indi[i], m, k, l);
            matrixSubtraction(block, 1, m, block_inv, 1, m, block);
            put_block_b(b, block, indi[i], m, k, l);
        }
    }

    matrixOutput(a, n, n, 3);

    for(i = 0; i < k; i++)
        for(j = 0; j < m; j++)
            x[indi[i] * m + j] = b[indi[i] * m + j];
    if(bl == k + 1)
        for(j = 0; j < l; j ++)
            x[indi[bl] * m + j] = b[indi[bl] * m + j];

    return 0;
}