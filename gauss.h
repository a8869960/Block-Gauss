//
// Created by varsem on 09.10.23.
//
#include <iostream>

using namespace std;

int inverseMatrix(double *a, double *A, double *B, int n, int *indi_m, int *indj_m, double norm); //a^(-1) = A

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
              double *block_h,
              double NORM);
double matrixNorm(double *A, int n);

void get_block(
        double* A,
        double* block,
        int i,
        int j,
        int n,
        int m,
        int k,
        int l);
void put_block(
        double* A,
        double* block,
        int i,
        int j,
        int n,
        int m,
        int k,
        int l);

void E(double* block, int m);

void get_block_b( double *B, double *block, int i, int m, int k, int l);
void put_block_b( double *B, double *block, int i, int m, int k, int l);

void matrix_product(double *A, double* B, double* C, int n, int s, int m);
