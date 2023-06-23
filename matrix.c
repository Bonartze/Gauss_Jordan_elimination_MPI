

#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

void fill_matrix(double *mtr, int N, int M) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++) {
            mtr[i * M + j] = rand() % 40;
            mtr[i * M + j]  = mtr[i * M + j] >-1 && mtr[i * M + j] <1 ? 3.3 : mtr[i * M + j];
        }
}

void printf_matrix(double *mtr, int N, int M) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++)
            printf("%.2f ", mtr[i * M + j]);
        printf("\n");
    }
}

void cat_matrices(double *o_matr, double *temp_matrix, double *matrix, int N, int M) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++)
            matrix[i * (N + M) + j] = temp_matrix[i * M + j];
        for (int j = M, t = 0; j < N + M; j++, t++) {
            matrix[i * (M + N) + j] = o_matr[i * N + t];
        }
    }
}