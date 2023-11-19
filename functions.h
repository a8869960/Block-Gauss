#include <iostream>

using namespace std;

int toInt(const char* str, int* ptr);

int fileMatrixInput(double* A, char* filename, int n);
void f(double *A, int s, int n);
void init_B(double *B, double *A, int n);

int gauss_func(int n,
               int m,
               double *A,
               [[maybe_unused]]double *B,
               double *x,
               int *indi_m,
               int *indj_m,
               int *indi,
               int *indj,
               double* a,
               double *b,
               double *block,
               double *block_inv,
               double *block_h);

void matrixOutput(double *matrix, int l, int n, int r);

int calc_r1(double* A, double* x, double* B, int n, double* helper, double *r1);
void calc_r2(double *x, int n, double *r2);

int matrixSubtraction(double* A1, double *A2, double *C, int n, int m);