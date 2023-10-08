#include <iostream>

using namespace std;

int toInt(const char* str, int* ptr);

int fileMatrixInput(double* A, char* filename, int n);
void f(double *A, int s, int n);
void init_B(double *B, double *A, int n);

int gauss_func(int n, int m, double *A, double *B, double *x, double *helper);

void matrixOutput(double *matrix, int l, int n, int r);

int calc_r1(double* A, double* x, double* B, int n, double* helper, double *r1);
void calc_r2(double *x, int n, double *r2);

int matrixProduct(double* A1,
                  int n,
                  int m,
                  double *A2,
                  int r,
                  int s,
                  double *C);
int matrixSubtraction(double* A1,
                      int n,
                      int m,
                      double *A2,
                      int r,
                      int s,
                      double *C);
void unit(double *a, double *b, double *c, int n, int m);

int inverseMatrix(double *a, double *A, double *B, int n, int *indi, int *indj);