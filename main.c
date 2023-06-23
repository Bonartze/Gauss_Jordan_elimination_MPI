#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include "matrix.h"

#define M 6
#define N 6

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int num_tasks;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    const int n_rows = N / num_tasks;

    int task_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    const int start_row = task_id * n_rows;
    const int end_row = start_row + n_rows;

    double temp_matrix[N * M];


    double m_chunk[(N + M) * n_rows];

    double pivot_row[N + M];
    double matrix[N * (N + M)];
    double o_matr[N * N];
    if (task_id == 0) {
        for (int i = 0; i < N * N; i++)
            o_matr[i] = 0.0;
        for (int i = 0; i < N; i++) {
            o_matr[N * i + i] = 1.0;
        }
        fill_matrix(temp_matrix, N, M);
        cat_matrices(o_matr, temp_matrix, matrix, N, M);
        printf("\n");
        printf_matrix(matrix, N, M + N);
    }
    MPI_Scatter(matrix, (N + M) * n_rows, MPI_DOUBLE, m_chunk,
                (N + M) * n_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Request requests[num_tasks];
    fill_matrix(temp_matrix, N, M);
    for (int row = 0; row < end_row; row++) {
        int mapped_rank = row / n_rows;
        if (task_id == mapped_rank) {
            // Calculate the row in the local matrix
            int local_row = row % n_rows;

            // Get the value of the pivot
            double pivot = m_chunk[local_row * (N + M) + row];

            // Divide the rest of the row by the pivot
            for (int col = row; col < (N + M); col++) {
                m_chunk[local_row * (N + M) + col] /= pivot;
            }

            // Send the pivot row to the other processes
            for (int i = mapped_rank + 1; i < num_tasks; i++) {
                MPI_Isend(m_chunk + (N + M) * local_row, M + N, MPI_DOUBLE, i, 0,
                          MPI_COMM_WORLD, &requests[i]);
            }

            // Eliminate the for the local rows
            for (int elim_row = local_row + 1; elim_row < n_rows; elim_row++) {
                // Get the scaling factor for elimination
                double scale = m_chunk[elim_row * (N + M) + row];

                // Remove the pivot
                for (int col = row; col < (N + M); col++) {
                    m_chunk[elim_row * (N + M) + col] -=
                            m_chunk[local_row * (N + M) + col] * scale;
                }
            }

            // Check if there are any outstanding messages
            for (int i = mapped_rank + 1; i < num_tasks; i++) {
                MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
            }
        } else {
            // Receive pivot row
            MPI_Recv(pivot_row, (N + M), MPI_DOUBLE, mapped_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            // Skip rows that have been fully processed
            for (int elim_row = 0; elim_row < n_rows; elim_row++) {
                // Get the scaling factor for elimination
                double scale = m_chunk[elim_row * (N + M) + row];

                // Remove the pivot
                for (int col = row; col < (N + M); col++) {
                    m_chunk[elim_row * (N + M) + col] -= pivot_row[col] * scale;
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int row = N - 1; row >= n_rows * task_id; row--) {
        double p_r[N + M];
        int mapped_rank = row / n_rows;
        if (task_id == mapped_rank) {
            int local_row = row % n_rows;
            for (int i = task_id - 1; i >= 0; i--) {
                MPI_Send(m_chunk + (N + M) * local_row, M + N, MPI_DOUBLE, i, 0,
                         MPI_COMM_WORLD);
            }
            for (int elim_row = local_row - 1; elim_row >= 0; elim_row--) {
                // Get the scaling factor for elimination
                double scale = m_chunk[elim_row * (M + N) + row];
                // Remove the pivot
                for (int col = row; col < M + N; col++) {
                    m_chunk[elim_row * (N + M) + col] -=
                            m_chunk[local_row * (N + M) + col] * scale;
                }
            }
        } else {
            MPI_Recv(p_r, (N + M), MPI_DOUBLE, mapped_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            for (int elim_row = 0; elim_row < n_rows; elim_row++) {
                double scale = m_chunk[elim_row * (N + M) + row];
                for (int col = row; col < (N + M); col++) {
                    m_chunk[elim_row * (N + M) + col] -= p_r[col] * scale;
                }
            }
        }
    }
    MPI_Gather(m_chunk, n_rows * (N + M), MPI_DOUBLE, matrix, n_rows * (N + M),
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (task_id == 0) {
        printf("\n");
        printf_matrix(matrix, N, (N + M));
    }
    MPI_Finalize();
    return 0;
}