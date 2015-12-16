#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>

#define N 4096 /* число точек */
#define h (1.0 / N) /* шаг сетки */
#define R 12 /* число редукций метода прогонки */
#define epsilon 1E-6
#define REPEATS 100

double my_func(double* y, double b, int n)
{
    if (n == 0) return y[0] - 1.0;
    else if (n == N) return y[N] - b;
    else
    {
        return (y[n + 1] - 2 * y[n] + y[n - 1]) / (h * h) - exp(y[n + 1]) / 12 - 5 * exp(y[n]) / 6 - exp(y[n - 1]) / 12;
    }
}

double norm(double* x)
{
    double ans = 0.0;
    int i;
    for(i = 0; i <= N; i++) ans += abs(x[i]);
    return ans;
}

double run(int threads_num, double b)
{
    int i, j, step,k;
    double* y = calloc(N + 1, sizeof(double)); /* сеточное решение */
    double* dy = calloc(N + 1, sizeof(double)); /* разность y^n-y^n+1 двух соседних приближений по итерациям метода Ньютона */
    double *A[R], *B[R], *C[R], *G[R]; /* коэффициенты трёхдиагональной системы для каждого шага редукции */
    double begin, end;
    omp_set_dynamic(0); /* нельзя динамически изменять количество нитей */
    omp_set_num_threads(threads_num); /* 4 нити */
    for(i = 0; i < R; i++)
    {
        A[i] = calloc(N + 1, sizeof(double));
        B[i] = calloc(N + 1, sizeof(double));
        C[i] = calloc(N + 1, sizeof(double));
        G[i] = calloc(N + 1, sizeof(double));
    }
    begin = omp_get_wtime(); /* начальная точка отсчёта времени */
    for( k = 0; k < REPEATS; k++){
        #pragma omp parallel private(i, j)
        {
                #pragma omp for
                for(i = 0; i <= N; i++) y[i] = 1.0 + (b - 1.0) * i / N; /* нулевое приближение */
                #pragma omp single
                {
                    dy[0] = dy[N] = 0.0;
                    for(j = 0; j < R; j++) B[j][0] = B[j][N] = 1.0; /* при редукции крайние значения матрицы одни и те же во всех итерациях метода Ньютона */
                }
                while(1) /* итерации метода Ньютона в цикле */
                {
                    #pragma omp for
                    for(i = 1; i < N; i++) /* изначальные значения коэффициентов */
                    {
                        B[0][i] = (-2.0 / (h * h) - 5 * exp(y[i]) / 6);
                        A[0][i] = (1.0 / (h * h) - exp(y[i - 1]) / 12);
                        C[0][i] = (1.0 / (h * h) - exp(y[i + 1]) / 12);
                        G[0][i] = my_func(y, b, i);
                    }
                    for(j = 1; j < R; j++) /* значения коэффициентов после редукции */
                    {
                        step = pow(2, j); /* шаг прогонки при редукции */
                        #pragma omp for
                        for(i = step; i < N; i += step)
                        {
                            B[j][i] = B[j - 1][i] - A[j - 1][i] * C[j - 1][i - step / 2] / B[j - 1][i - step / 2] - C[j - 1][i] * A[j - 1][i + step / 2] / B[j - 1][i + step / 2];
                            A[j][i] = - A[j - 1][i] * A[j - 1][i - step / 2] / B[j - 1][i - step / 2];
                            C[j][i] = - C[j - 1][i] * C[j - 1][i + step / 2] / B[j - 1][i + step / 2];
                            G[j][i] = G[j - 1][i] - A[j - 1][i] * G[j - 1][i - step / 2] / B[j - 1][i - step / 2] - C[j - 1][i] * G[j - 1][i + step / 2] / B[j - 1][i + step / 2];
                        }
                    } /* редукция прогонки завершена */
                    #pragma omp single
                    {
                        dy[N / 2] = G[R - 1][N / 2] / B[R - 1][N / 2]; /* первый обратный шаг редукции */
                        dy[N / 4] = (G[R - 2][N / 4] - C[R - 2][N / 4] * dy[N / 2]) / B[R - 2][N / 4];
                        dy[N * 3 / 4] = (G[R - 2][N * 3 / 4] - A[R - 2][N * 3 / 4] * dy[N / 2] ) / B[R - 2][N * 3 / 4]; /* второй обратный шаг редукции */
                    }
                    for(j = R - 3; j >= 0; j--)
                    {
                        step = pow(2, j);
                        #pragma omp for
                        for(i = step; i < N; i += 2 * step) dy[i] = (G[j][i] - C[j][i] * dy[i + step] - A[j][i] * dy[i - step]) / B[j][i];
                    } /* оставшиеся обратные шаги редукции */
                    #pragma omp for
                    for(i = 0; i <= N; i++) y[i] -= dy[i]; /* одна итерация метода Ньютона */
                    if (norm(dy) < epsilon) break; /* условие останова метода Ньютона */
                }
        }
    }
    end = omp_get_wtime(); /* конечная точка отсчёта времени */
                    for(i = 0; i < R; i++)
                    {
                        free(A[i]);
                        free(B[i]);
                        free(C[i]);
                        free(G[i]);
                    }
    if( threads_num == 1){
        char str_dest[50];
        sprintf( str_dest, "prog_1_b_%f_results.txt",b);
        FILE* fp = fopen(str_dest, "w"); /* вывод полученной функции в файл */
        fprintf(fp, "X\tY\r\n");
        for(i = 0; i <= N; i++) fprintf(fp, "%e\t%e\r\n", ((double) i / N), y[i]);
        fclose(fp);
    }
    free(y);
    free(dy);
    return (end - begin)/REPEATS;
}





int main(int argc, char** argv)
{
    double b = atof(argv[1]);
    double T1 = run(1, b); /* время выполнения одного треда*/
    char str_dest[50];
    sprintf( str_dest, "prog_1_b_%f_time.txt",b);
    FILE* fp = fopen( str_dest, "w");
    fprintf(fp, "N S E\r\n");
    int i;
    double Ti;
    for(i = 2; i <= 16; i++)
    {
        Ti = run(i, b);
        fprintf(fp, "%d %f %f\r\n", i, (T1 / Ti), (T1 / Ti / i));
    }
    fclose(fp);
    return 0;
}

