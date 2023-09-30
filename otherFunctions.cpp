//
// Created by varsem on 28.09.23.
//
#define x(i) x[i - 1]
#define helper(i) helper[i - 1]

#include <stdlib.h>
#include <errno.h>
#include <limits.h>

#include <cstring>

#include "functions.h"

int toInt(const char* str, int* ptr)
{
    long L;
    char* e;

    errno = 0;
    L = strtol(str, &e, 10);

    if (!errno && *e == '\0')
        if (INT_MIN <= L && L <= INT_MAX)
        {
            *ptr = (int)L;
            return 0;
        }
        else
            return -1;
    else
        return -1;
}

double norm1(double *x, int n)
{
    double norma = 0;

    for(int i = 1; i <= n; i++)
        norma += abs(x(i));

    return norma;
}

int calc_r1(double* A, double* x, double* B, int n, double* helper, double *r1)
{
    memset(helper, 0, n);

    if(matrixProduct(A, n, n, x, n, 1, helper) == -1)
        return -1;

    if(matrixSubtraction(helper, n, 1, B, n, 1, helper) == -1)
        return -1;

    *r1 = norm1(helper, n) / norm1(B, n);

    return 0;
}

void calc_r2(double *x, int n, double *r2)
{
    for(int i = 1; i <= n; i++)
        *r2 += abs(x(i) - (i % 2));
}