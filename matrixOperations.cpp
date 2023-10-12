//
// Created by varsem on 04.10.23.
//
#include <iostream>

using namespace std;

int matrixProduct(double* A1,
                  int n,
                  int m,
                  double *A2,
                  int r,
                  int s,
                  double *C)
{
    if(m != r)
    {
        cout << "Can't calculate the product of matrix." << endl;
        return -1;
    }

    for(int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= s; j++)
        {
            C[s * (i - 1) + j - 1] = 0;
            for(int k = 1; k <= m; k++)
            {
                C[s * (i - 1) + j - 1] += A1[m * (i - 1) + k - 1] * A2[s * (k - 1) + j - 1];
            }
        }

    }
    return 0;
}

int matrixSubtraction(double* A1,
                      int n,
                      int m,
                      double *A2,
                      int r,
                      int s,
                      double *C)
{
    if(n != r or m != s)
    {
        cout << "Can't do matrix subtraction." << endl;
        return -1;
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
            C[m * i + j] = A1[m * i + j] - A2[m * i + j];
    }

    return 0;
}

void unit(double *a, double *b, double *c, int n, int m)
{
    int k, l, bl, i, j, v, h, r, t, ah, v3, h3, q, s;
    double *pa, *pb, *pc;
    double s00, s01, s02, s10, s11, s12, s20, s21, s22,  sum;

    k = n / m;
    l = n - k * m;
    bl = (l != 0) ? k + 1 : k;

    for(i = 0; i < bl; i++)
        for(j = 0; j < bl; j++)
        {
            v = (i < k) ? m : l;
            h = (j < k) ? m : l;
            pc = c + i * n * m + j * m;

            for(r = 0; r < v; r++)
                for(t = 0; t < h; t++)
                pc[r * n + t] = 0;

            for(s = 0; s < bl; s++)
            {
                ah = (s < k) ? m : l;
                pa = a + i * n * m + s * m;
                pb = b + s * n * m + j * m;

                v3 = v % 3;
                h3 = h % 3;

                for(r = 0; r < v3; r++)
                {
                    for (t = 0; t < h3; t++)
                    {
                        sum = 0;

                        for(q = 0; q < ah; q++)
                            sum += pa[r * n + q] * pb[q * n + t];
                        pc[r * n + t] += s;
                    }

                    for(; t < h; t += 3)
                    {
                        s00 = 0; s01 = 0; s02 = 0;

                        for(q = 0; q < ah; q++)
                        {
                            s00 += pa[r * n + q] * pb[q * n + t];
                            s01 += pa[r * n + q] * pb[q * n + t + 1];
                            s02 += pa[r * n + q] * pb[q * n + t + 2];
                        }

                        pc[r * n + t] += s00;
                        pc[r * n + t + 1] += s01;
                        pc[r * n + t + 2] += s02;
                    }
                }

                for(; r < v; r += 3)
                {
                    for(t = 0; t < h3; t++)
                    {
                        s00 = 0;
                        s10 = 0;
                        s20 = 0;

                        for(q = 0; q < h; q++)
                        {
                            s00 += pa[r * n + q] * pb[q * n + t];
                            s10 += pa[(r + 1) * n + q] * pb[q * n + t];
                            s20 += pa[(r + 2) * n + q] * pb[q * n + t];
                        }

                        pc[r * n + t] += s00;
                        pc[(r + 1) * n + t + 1] += s10;
                        pc[(r + 2) * n + t + 2] += s20;
                    }

                    for(; t < h; t += 3)
                    {
                        s00 = 0; s01 = 0; s02 = 0;
                        s10 = 0; s11 = 0; s12 = 0;
                        s20 = 0; s21 = 0; s22 = 0;

                        for(q = 0; q < ah; q++)
                        {
                            s00 += pa[r * n + q] * pb[q * n + t];
                            s01 += pa[r * n + q] * pb[q * n + t + 1];
                            s02 += pa[r * n + q] * pb[q * n + t + 2];
                            s10 += pa[(r + 1) * n + q] * pb[q * n + t];
                            s11 += pa[(r + 1) * n + q] * pb[q * n + t + 1];
                            s12 += pa[(r + 1) * n + q] * pb[q * n + t + 2];
                            s20 += pa[(r + 2) * n + q] * pb[q * n + t];
                            s21 += pa[(r + 2) * n + q] * pb[q * n + t + 1];
                            s22 += pa[(r + 2) * n + q] * pb[q * n + t + 2];
                        }

                        pc[r * n + t] += s00;
                        pc[r * n + t + 1] += s01;
                        pc[r * n + t + 2] += s02;
                        pc[(r + 1) * n + t] += s10;
                        pc[(r + 1) * n + t + 1] += s11;
                        pc[(r + 1) * n + t + 2] += s12;
                        pc[(r + 2) * n + t] += s20;
                        pc[(r + 2) * n + t + 1] += s21;
                        pc[(r + 2) * n + t + 2] += s22;
                    }
                }
            }
        }
}

void matrix_product(double *A, double* B, double* C, int n, int s, int m)
{
    int k1, l1, k2, l2, k3, l3, i, j, t;
    double t00, t01, t02, t10, t11, t12, t20, t21, t22;
    k1 = n/3;
    k2 = s/3;
    k3 = m/3;
    l1 = n - k1*3;
    l2 = s - k2*3;
    l3 = m - k3*3;

    for(i = 0; i < k1; i++){
        for(j = 0; j < k3; j++){
            t00 = 0;
            t01 = 0;
            t02 = 0;
            t10 = 0;
            t11 = 0;
            t12 = 0;
            t20 = 0;
            t21 = 0;
            t22 = 0;
            for(t = 0; t < k2; t++){
                t00 += A[3*i*s+3*t]*B[3*t*m+3*j] + A[3*i*s+3*t+1]*B[(3*t+1)*m+3*j] + A[3*i*s+3*t+2]*B[(3*t+2)*m+3*j];
                t01 += A[3*i*s+3*t]*B[3*t*m+3*j+1] + A[3*i*s+3*t+1]*B[(3*t+1)*m+3*j+1] + A[3*i*s+3*t+2]*B[(3*t+2)*m+3*j+1];
                t02 += A[3*i*s+3*t]*B[3*t*m+3*j+2] + A[3*i*s+3*t+1]*B[(3*t+1)*m+3*j+2] + A[3*i*s+3*t+2]*B[(3*t+2)*m+3*j+2];
                t10 += A[(3*i+1)*s+3*t]*B[3*t*m+3*j] + A[(3*i+1)*s+3*t+1]*B[(3*t+1)*m+3*j] + A[(3*i+1)*s+3*t+2]*B[(3*t+2)*m+3*j];
                t11 += A[(3*i+1)*s+3*t]*B[3*t*m+3*j+1] + A[(3*i+1)*s+3*t+1]*B[(3*t+1)*m+3*j+1] + A[(3*i+1)*s+3*t+2]*B[(3*t+2)*m+3*j+1];
                t12 += A[(3*i+1)*s+3*t]*B[3*t*m+3*j+2] + A[(3*i+1)*s+3*t+1]*B[(3*t+1)*m+3*j+2] + A[(3*i+1)*s+3*t+2]*B[(3*t+2)*m+3*j+2];
                t20 += A[(3*i+2)*s+3*t]*B[3*t*m+3*j] + A[(3*i+2)*s+3*t+1]*B[(3*t+1)*m+3*j] + A[(3*i+2)*s+3*t+2]*B[(3*t+2)*m+3*j];
                t21 += A[(3*i+2)*s+3*t]*B[3*t*m+3*j+1] + A[(3*i+2)*s+3*t+1]*B[(3*t+1)*m+3*j+1] + A[(3*i+2)*s+3*t+2]*B[(3*t+2)*m+3*j+1];
                t22 += A[(3*i+2)*s+3*t]*B[3*t*m+3*j+2] + A[(3*i+2)*s+3*t+1]*B[(3*t+1)*m+3*j+2] + A[(3*i+2)*s+3*t+2]*B[(3*t+2)*m+3*j+2];
            }
            for(t = 0; t < l2; t++){
                t00 += A[3*i*s+t+3*k2]*B[(t+3*k2)*m+3*j];
                t01 += A[3*i*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+1];
                t02 += A[3*i*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+2];
                t10 += A[(3*i+1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j];
                t11 += A[(3*i+1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+1];
                t12 += A[(3*i+1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+2];
                t20 += A[(3*i+2)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j];
                t21 += A[(3*i+2)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+1];
                t22 += A[(3*i+2)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+2];
            }
            C[3*i*m+3*j] = t00;
            C[3*i*m+3*j+1] = t01;
            C[3*i*m+3*j+2] = t02;
            C[(3*i+1)*m+3*j] = t10;
            C[(3*i+1)*m+3*j+1] = t11;
            C[(3*i+1)*m+3*j+2] = t12;
            C[(3*i+2)*m+3*j] = t20;
            C[(3*i+2)*m+3*j+1] = t21;
            C[(3*i+2)*m+3*j+2] = t22;
        }

        for(j = 0; j < l3; j++){
            t00 = 0;
            t10 = 0;
            t20 = 0;
            for(t = 0; t < k2; t++){

                t00 += A[3*i*s+3*t]*B[3*t*m+j+3*k3] + A[3*i*s+3*t+1]*B[(3*t+1)*m+j+3*k3] + A[3*i*s+3*t+2]*B[(3*t+2)*m+j+3*k3];
                t10 += A[(3*i+1)*s+3*t]*B[3*t*m+j+3*k3] + A[(3*i+1)*s+3*t+1]*B[(3*t+1)*m+j+3*k3] + A[(3*i+1)*s+3*t+2]*B[(3*t+2)*m+j+3*k3];
                t20 += A[(3*i+2)*s+3*t]*B[3*t*m+j+3*k3] + A[(3*i+2)*s+3*t+1]*B[(3*t+1)*m+j+3*k3] + A[(3*i+2)*s+3*t+2]*B[(3*t+2)*m+j+3*k3];

            }
            for(t = 0; t < l2; t++){
                t00 += A[3*i*s+t+3*k2]*B[(t+3*k2)*m+j+3*k3];
                t10 += A[(3*i+1)*s+(t+3*k2)]*B[(t+3*k2)*m+j+3*k3];
                t20 += A[(3*i+2)*s+(t+3*k2)]*B[(t+3*k2)*m+j+3*k3];
            }
            C[3*i*m+j+3*k3] = t00;
            C[(3*i+1)*m+j+3*k3] = t10;
            C[(3*i+2)*m+j+3*k3] = t20;
        }
    }

    for(i = 0; i < l1; i++){
        for(j = 0; j < k3; j++){
            t00 = 0;
            t01 = 0;
            t02 = 0;
            for(t = 0; t < k2; t++){
                t00 += A[(i+3*k1)*s+3*t]*B[3*t*m+3*j] + A[(i+3*k1)*s+3*t+1]*B[(3*t+1)*m+3*j] + A[(i+3*k1)*s+3*t+2]*B[(3*t+2)*m+3*j];
                t01 += A[(i+3*k1)*s+3*t]*B[3*t*m+3*j+1] + A[(i+3*k1)*s+3*t+1]*B[(3*t+1)*m+3*j+1] + A[(i+3*k1)*s+3*t+2]*B[(3*t+2)*m+3*j+1];
                t02 += A[(i+3*k1)*s+3*t]*B[3*t*m+3*j+2] + A[(i+3*k1)*s+3*t+1]*B[(3*t+1)*m+3*j+2] + A[(i+3*k1)*s+3*t+2]*B[(3*t+2)*m+3*j+2];

            }
            for(t = 0; t < l2; t++){
                t00 += A[(i+3*k1)*s+t+3*k2]*B[(t+3*k2)*m+3*j];
                t01 += A[(i+3*k1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+1];
                t02 += A[(i+3*k1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+2];
            }
            C[(i+3*k1)*m+3*j] = t00;
            C[(i+3*k1)*m+3*j+1] = t01;
            C[(i+3*k1)*m+3*j+2] = t02;
        }
        for(j = 0; j < l3; j++){
            t00 = 0;
            for(t = 0; t < k2; t++){
                t00 += A[(i+3*k1)*s+3*t]*B[3*t*m+j+3*k3] + A[(i+3*k1)*s+3*t+1]*B[(3*t+1)*m+j+3*k3] + A[(i+3*k1)*s+3*t+2]*B[(3*t+2)*m+j+3*k3];
            }
            for(t = 0; t < l2; t++){
                t00 += A[(i+3*k1)*s+t+3*k2]*B[(t+3*k2)*m+j+3*k3];
            }
            C[(i+3*k1)*m+j+3*k3] = t00;
        }
    }
}
