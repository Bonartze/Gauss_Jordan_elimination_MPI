#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include "matrix.h"

#define M 4
#define N 3

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int num_tasks;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    // Calculate the number of rows mapped to each process
    // Assumes this divides evenly
    const int n_rows = N / num_tasks;

    // Get the task ID
    int task_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    const int start_row = task_id * n_rows;
    const int end_row = start_row + n_rows;

    // Matrix - Only initialized in rank 0
    double temp_matrix[N * M];


    // Each process will store a chunk of the matrix
    double m_chunk[(N + M) * n_rows];

    // Each process will receive a pivot row each iteration
    double pivot_row[N + M];
    double matrix[N * (N + M)];
    // Only rank 0 create initializes the matrix
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
    // Before doing anything, send parts of the matrix to each process
    MPI_Scatter(matrix, (N + M) * n_rows, MPI_DOUBLE, m_chunk,
                (N + M) * n_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Store requests that for non-blocking sends
    MPI_Request requests[num_tasks];

    // Performance gaussian elimination
    for (int row = 0; row < end_row; row++) {
        // See if this process is responsible for the pivot calculation
        int mapped_rank = row / n_rows;

        // If the row is mapped to this rank
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
    MPI_Request reqs[num_tasks];
    int r_end_row = task_id - 1;
    for (int row = N - 1; row > r_end_row; row--) {
        double p_r[N + M];
        printf("%i: >%i<\n", task_id, row);
        int mapped_rank = row / n_rows;
        if (task_id == mapped_rank) {
            int local_row = row % n_rows;
            //  double pivot = m_chunk[local_row * (N + M) + row + 1];
            // if (row != N - 1) {
            //     for (int col = 0; col < N + M; col++) {
            //         m_chunk[local_row * (N + M) + col] /= pivot;
            //     }
            // }
            /*if (task_id == 2)
                for (int i = 0; i < M + N; i++)
                    printf("%f ", m_chunk[i]);*/
            for (int i = task_id - 1; i >= 0; i--) {
                MPI_Send(m_chunk + (N + M) * local_row, M + N, MPI_DOUBLE, i, 0,
                         MPI_COMM_WORLD);
            }
            /*if (row != N - 1) {
                for (int elim_row = local_row - 1; elim_row > 0; elim_row--) {
                    // Get the scaling factor for elimination
                    double scale = m_chunk[elim_row * (N + M) + row];

                    // Remove the pivot
                    for (int col = row; col < N + M; col++) {
                        m_chunk[elim_row * (N + M) + col] -=
                                m_chunk[local_row * (N + M) + col] * scale;
                    }
                }
            }*/
            /*for (int i = mapped_rank - 1; i > -1; i--) {
                MPI_Wait(&reqs[i], MPI_STATUS_IGNORE);
            }*/
        } else {/*
            if (task_id == 0) {
                printf(">>%i<<\n", row);
            }*/
            /*for (int i = task_id + 1; i < num_tasks; i++)
                MPI_Wait(&reqs[i], MPI_STATUS_IGNORE);*/
            MPI_Recv(p_r, (N + M), MPI_DOUBLE, mapped_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            for (int elim_row = 0; elim_row < n_rows; elim_row++) {
                // Get the scaling factor for elimination
                double scale = m_chunk[elim_row * (N + M) + row];
                // Remove the pivot
                for (int col = row; col < (N + M); col++) {
                    m_chunk[elim_row * (N + M) + col] -= p_r[col] * scale;
                }
            }
        }
    }
    // printf("%i: Almost...\n", task_id);
    MPI_Gather(m_chunk, n_rows * (N + M), MPI_DOUBLE, matrix, n_rows * (N + M),
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (task_id == 0) {
        printf("\n");
        printf_matrix(matrix, N, (N + M));
    }
    MPI_Finalize();
    return 0;
}