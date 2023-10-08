//
// Created by varsem on 04.10.23.
//
#define eps 1e-20
#include <iostream>
#include <cstring>

#include "functions.h"

using namespace std;

//int inverseMatrix(double *a, double *A, double *B, int n, int *indi, int *indj);
int matrixMax(double *A, int k, int n, int *indi, int *indj);

int inverseMatrix(double *a, double *A, double *B, int n, int *indi, int *indj)
{
    int i, j, ii;

    for(i = 0; i < n; i++)
    {
        indi[i] = i;
        indj[i] = i;
    }

    memcpy(A, a, sizeof(double) * n * n);

//    cout << "--------Ah---" << endl;
//    matrixOutput(A, n, n, n);

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
        {
            if(i == j)
                B[i * n + j] = 1;
            else
                B[i * n + j] = 0;
        }

    //Прямой ход
    for(i = 0; i < n; i++)
    {
        if(matrixMax(A, i, n, indi, indj) == -1)
            return -1;

        double Aii = A[indi[i] * n + indj[i]];
//        cout << "A" << i << " " << Aii << endl;

        for(j = 0; j < i + 1; j++)
            B[indi[i] * n + indj[j]] = B[indi[i] * n + indj[j]] / Aii;
        for(j = i + 1; j < n; j++)
        {
            A[indi[i] * n + indj[j]] = A[indi[i] * n + indj[j]] / Aii;
            B[indi[i] * n + indj[j]] = B[indi[i] * n + indj[j]] / Aii;
        }
        A[indi[i] * n + indj[i]] = 1;

//        cout << "--------A/---" << endl;
//        matrixOutput(A, n, n, n);

        for(ii = i + 1; ii < n; ii++)
        {
            for(j = 0; j < i + 1; j++)
                B[indi[ii] * n + indj[j]] = B[indi[ii] * n + indj[j]] - B[indi[i] * n + indj[j]] * A[indi[ii] * n + indj[i]];
            for(j = i + 1; j < n; j++)
            {
                A[indi[ii] * n + indj[j]] = A[indi[ii] * n + indj[j]] - A[indi[i] * n + indj[j]] * A[indi[ii] * n + indj[i]];
                B[indi[ii] * n + indj[j]] = B[indi[ii] * n + indj[j]] - B[indi[i] * n + indj[j]] * A[indi[ii] * n + indj[i]];
            }
//            B[indi[ii] * n + indj[i]] = B[indi[ii] * n + indj[i]] - B[indi[i] * n + indj[i]] * A[indi[ii] * n + indj[i]];
            A[indi[ii] * n + indj[i]] = 0;
        }
//        cout << "--------A0---" << endl;
//        matrixOutput(A, n, n, n);
    }
//    cout << "Pr hod A" << endl;
//    matrixOutput(A, n, n, n);
//    cout << "Pr hod B" << endl;
//    matrixOutput(B, n, n, n);

    //Обратный ход
    for(i = n - 1; i >= 0; i--)
    {
        for(ii = i - 1; ii >= 0; ii--)
        {
            double Aa = A[indi[ii] * n + indj[i]];

            for(j = n - 1; j >= 0; j--)
            {
                B[indi[ii] * n + indj[j]] = B[indi[ii] * n + indj[j]] - Aa * B[indi[i] * n + indj[j]];
                A[indi[ii] * n + indj[j]] = A[indi[ii] * n + indj[j]] - Aa * A[indi[i] * n + indj[j]];
            }
        }
    }

//    memcpy(A, B, sizeof(double) * n * n);
//    for(i = 0; i < n; i++)
//    {
//        for(j = 0; j < n; j++)
//        {
//            cout << " " << A[indi[i] * n + indj[j]];
//        }
//        cout << endl;
//    }

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            A[indi[i] * n + j] = B[indj[i] * n + j];

//    cout << "Ob hod B" << endl;
//    matrixOutput(B, n, n, n);
//
    cout << "Ob hod A" << endl;
    matrixOutput(A, n, n, n);

//    cout << "Matrix product" << endl;
//    unit(B, a, A, n, n);
//    matrixProduct(A, n, n, a, n, n, B);
//    matrixOutput(B, n, n, n);

    return 0;
}

int matrixMax(double *A, int k, int n, int *indi, int *indj)
{
    double max = 0;
    int imax = k, jmax = k;

    for(int i = k; i < n; i++)
        for(int j = k; j < n; j++)
        {
            if(abs(A[indi[i] * n + indj[j]]) > max)
            {
                max = abs(A[indi[i] * n + indj[j]]);
                imax = i;
                jmax = j;
            }
        }

//    cout << "MAX = " << imax << jmax << " " << A[indi[imax] * n + indj[jmax]] << endl;

    if(abs(max) < eps)
        return -1;

    int helper;

    helper = indi[imax];
    indi[imax] = indi[k];
    indi[k] = helper;

    helper = indj[jmax];
    indj[jmax] = indj[k];
    indj[k] = helper;

    return 0;
}