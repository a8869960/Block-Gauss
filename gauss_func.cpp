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

//        cout << "BLOCK II" << endl;
//        matrixOutput(block_inv, m, m, 3);

        //Деление первой строчки
        E(block, m);
        put_block(a, block, indi[step], indj[step], n, m, k, l);
        for(j = step + 1; j < bl; j++)
        {
            get_block(a, block, indi[step], indj[j], n, m, k, l);
//            cout << "BLOCK" << endl;
//            matrixOutput(block, m, m, 3);
            matrix_product(block_inv, block, block_h, m, m, m);
//            cout << "BLOCK H" << endl;
//            matrixOutput(block_h, m, m, 3);
            put_block(a, block_h, indi[step], indj[j], n, m, k, l);
        }
        get_block_b(b, block, indi[step], m, k, l);
        matrix_product(block_inv, block, block_h, m, m, 1);
        put_block_b(b, block_h, indi[step], m, k, l);

//        cout << "A /" << endl;
//        matrixOutput(a, n, n, 3);
//        cout << "B /" << endl;
//        matrixOutput(b, 1, n, 3);
//        cout << endl;

        //Зануление столбца
        for(i = step + 1; i < bl; i++)
        {
            get_block(a, block, indi[i], indj[step], n, m, k ,l);
//            cout << "BLOCK" << endl;
//            matrixOutput(block, m, m, 3);
            for(j = step + 1; j < bl; j++)
            {
                get_block(a, block_h, indi[step], indj[j], n, m, k, l);

//                cout << "BLOCK_H" << endl;
//                matrixOutput(block_h, m, m, 3);

                matrix_product(block, block_h, block_inv, m, m, m);

//                cout << "BLOCK_INV" << endl;
//                matrixOutput(block_inv, m, m, 3);

                get_block(a, block_h, indi[i], indj[j], n, m, k, l);
//                cout << "BLOCK_H" << endl;
//                matrixOutput(block_h, m, m, 3);
                matrixSubtraction(block_h, m, m, block_inv, m, m, block_h);
//                cout << "BLOCK_H" << endl;
//                matrixOutput(block_h, m, m, 3);
                put_block(a, block_h, indi[i], indj[j], n, m, k, l);
            }
            get_block_b(b, block_h, indi[step], m, k, l);
//            cout << "BLOCK_H" << endl;
//            matrixOutput(block_h, m, m, 3);
            matrix_product(block, block_h, block_inv, m, m, 1);

//            matrixProduct(block, m, m, block_h, m, 1, block_inv);
            get_block_b(b, block_h, indi[i], m, k, l);
            matrixSubtraction(block_h, 1, m, block_inv, 1, m, block_h);
            put_block_b(b, block_h, indi[i], m, k, l);
        }
        memset(block, 0, sizeof(double) * m * m);
        for(i = step + 1; i < bl; i++)
            put_block(a, block, indi[i], indj[step], n, m, k, l);

//        cout << "A -" << endl;
//        matrixOutput(a, n, n, 3);
//        cout << "B -" << endl;
//        matrixOutput(b, 1, n, 3);
//        cout << endl;
    }

    if(bl == k + 1)
    {
        get_block(a, block, indi[k], indj[k], n, m, k, l);
        inverseMatrix(block, block_inv, block_h, l, indi_m, indj_m);

        E(block, l);
        put_block(a, block, indi[k], indj[k], n, m, k, l);
//        cout << "PR HOD A" << endl;
//        matrixOutput(a, n, n, 5);

        get_block_b(b, block, indi[k], m, k, l);
        matrix_product(block_inv, block, block_h, l, l, 1);
        put_block_b(b, block_h, indi[k], m, k, l);
    }

//    cout << "PR HOD A" << endl;
//    matrixOutput(a, n, n, 5);
//    cout << "PR HOD B" << endl;
//    matrixOutput(b, 1, n, 5);
//    cout << endl;

    //Обратный ход
    if(bl == k + 1)
    {
        get_block_b(b, block, indi[k], m, k, l);
//        cout << "BLOCK" << endl;
//        matrixOutput(block, m, m, 3);

        for(i = k - 1; i >= 0; i--)
        {
            get_block(a, block_h, indi[i], indj[k], n, m, k, l);
//            cout << "BLOCK_H" << endl;
//            matrixOutput(block_h, m, m, 3);

            matrix_product(block_h, block, block_inv, l, l, 1);
//            cout << "BLOCK_INV" << endl;
//            matrixOutput(block_inv, m, m, 3);

            get_block_b(b, block_h, indi[i], m, k, l);
//            cout << "BLOCK_H" << endl;
//            matrixOutput(block_h, 1, m, m);

            matrixSubtraction(block_h, 1, l, block_inv, 1, l, block_h);
//            cout << "BLOCK_H" << endl;
//            matrixOutput(block_h, 1, m, m);

            put_block_b(b, block_h, indi[i], m, k, l);
        }
    }
    for(step = k - 1; step >= 0; step--)
    {
        get_block_b(b, block, indi[step], m, k, l);

        for(i = step - 1; i >= 0; i--)
        {
            get_block(a, block_h, indi[i], indj[step], n, m, k, l);
            matrix_product(block_h, block, block_inv, m, m, 1);
            get_block_b(b, block_h, indi[i], m, k, l);
            matrixSubtraction(block_h, 1, m, block_inv, 1, m, block_h);
            put_block_b(b, block_h, indi[i], m, k, l);
        }
    }

//    cout << "OB HOD A" << endl;
//    matrixOutput(a, n, n, 5);
//    cout << "OB HOD B" << endl;
//    matrixOutput(b, 1, n, 5);
//    cout << endl;

    memcpy(x, b, sizeof(double) * n);

    return 0;
}