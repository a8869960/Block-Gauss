#include <ctime>

#include "functions.h"

int main(int ac, char* av[])
{
    try
    {
        int n, m, r, s;
        char* filename = nullptr;
        double *A = nullptr;

        if(ac != 5 and ac != 6)
        {
            cout << "Wrong parameters." << endl;
            return -1;
        }

        if(toInt(av[1], &n) == -1)
            return -1;

        if(toInt(av[2], &m) == -1)
            return -1;

        if(toInt(av[3], &r) == -1)
            return -1;

        if(toInt(av[4], &s) == -1)
            return -1;

        if(s == 0 and ac == 6)
        {
            filename = av[5];
        } else if((s == 0 and ac == 5) or (s != 0 and ac == 6))
        {
            cout << "Wrong parameter s and filename." << endl;
            return -1;
        }

        //Проверка аргументов командной строки
        if(n < 1 or m < 1 or n < m)
        {
            cout << "Wrong n and m parameters." << endl;
            return -1;
        }
        if(s < 0 or s > 4)
        {
            cout << "Wrong parameter s." << endl;
            return -1;
        }
        if(r < 1)
        {
            cout << "Wrong parameter r." << endl;
            return -1;
        }

        //Инициализация матрицы A
        A = new double[n*n];

        if(ac == 5)
        {
            f(A, s, n);
        } else
        {
            if(fileMatrixInput(A, filename, n) == -1)
            {
                delete[] A;
                return -1;
            }
        }

        cout << "Matrix A:" << endl;
        matrixOutput(A, n, n, r);

        double *B = new double[n];
        init_B(B, A, n);

        cout << "Matrix B:" << endl;
        matrixOutput(B, 1, n, r);

        double *x = new double[n];
        double *helper = nullptr;

        //Метод Гаусса
        clock_t start_time =  clock();
        if(gauss_func(n, m, A, B, x, helper) == -1)
        {
            cout << "Can't be used this method." << endl;
        } else
        {
            cout << "Solution:" << endl;
            matrixOutput(x, 1, n, r);
        }
        double t1 = (start_time - clock()) / CLOCKS_PER_SEC;

        //Запись результата
        start_time = clock();
        double r1 = 0;
        helper = new double[n];
        if(calc_r1(A, x, B, n, helper, &r1) == -1)
            return -1;

        double r2;
        calc_r2(x, n, &r2);
        double t2 = (start_time - clock()) / CLOCKS_PER_SEC;

        int task = 11;

        printf (
                "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
                av[0], task, r1, r2, t1, t2, s, n, m);

        delete[] A;
        delete[] B;
        delete[] x;
        delete[] helper;

        return 0;
    } catch (const bad_alloc& e)
    {
        cout << "Bad alloc" << endl;
        return -1;
    }
}