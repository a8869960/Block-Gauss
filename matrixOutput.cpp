//
// Created by varsem on 27.09.23.
//
#define matrix(i, j) matrix[(i - 1) * n + j - 1]

#include <cstdio>

#include "functions.h"

int min(int r, int l);

void matrixOutput(double *matrix, int l, int n, int r)
{

    for(int i = 1; i <= min(r, l); i++)
    {
        for(int j = 1; j <= r; j++)
        {
            printf(" %10.3e", matrix(i, j));
        }
        printf("\n");
    }
}

int min(int r, int l)
{
    return (r < l) ? r : l;
}