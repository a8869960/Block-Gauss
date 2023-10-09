//
// Created by varsem on 09.10.23.
//
#include <iostream>

int inverseMatrix(double *a, double *A, double *B, int n, int *indi_m, int *indj_m); //a^(-1) = A

int matrixMax(double *A, int n, int m, int k, int *indi, int *indj);
double matrixNorm(double *A, int n);

void get_block(
        double* A,
        double* block,
        size_t i,
        size_t j,
        size_t n,
        size_t m,
        size_t k,
        size_t l);
void put_block(
        double* A,
        double* block,
        size_t i,
        size_t j,
        size_t n,
        size_t m,
        size_t k,
        size_t l);
