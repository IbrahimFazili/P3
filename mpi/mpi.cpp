// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// // Global variables
// double *h_local = nullptr, *u_local = nullptr, *v_local = nullptr;
// int nx_local = 0, ny_local = 0;
// int rank = 0, num_procs = 0;
// double dt, g;  // Global gravity constant and time step

// // Initialize function
// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_) {
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

//     printf("Rank %d: Initializing with nx = %d, ny = %d\n", rank, nx_, ny_);
    
//     // Store global constants
//     dt = dt_;  // Initialize the global time step
//     g = g_;    // Initialize the global gravity constant

//     // Calculate local domain size
//     nx_local = nx_ / num_procs;
//     ny_local = ny_;

//     // Allocate local arrays and check allocation success
//     h_local = (double*)malloc(nx_local * ny_local * sizeof(double));
//     u_local = (double*)malloc(nx_local * ny_local * sizeof(double));
//     v_local = (double*)malloc(nx_local * ny_local * sizeof(double));
//     if (!h_local || !u_local || !v_local) {
//         fprintf(stderr, "Rank %d: Memory allocation failed\n", rank);
//         MPI_Abort(MPI_COMM_WORLD, -1);
//     }
//     printf("Rank %d: Memory allocated for h_local, u_local, v_local\n", rank);

//     // Scatter data if rank is 0 (only rank 0 has the full arrays h0, u0, v0)
//     if (rank == 0) {
//         MPI_Scatter(h0, nx_local * ny_local, MPI_DOUBLE, h_local, nx_local * ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatter(u0, nx_local * ny_local, MPI_DOUBLE, u_local, nx_local * ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatter(v0, nx_local * ny_local, MPI_DOUBLE, v_local, nx_local * ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         free(h0);
//         free(u0);
//         free(v0);
//         printf("Rank %d: Data scattered from root process\n", rank);
//     } else {
//         MPI_Scatter(nullptr, nx_local * ny_local, MPI_DOUBLE, h_local, nx_local * ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatter(nullptr, nx_local * ny_local, MPI_DOUBLE, u_local, nx_local * ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatter(nullptr, nx_local * ny_local, MPI_DOUBLE, v_local, nx_local * ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         printf("Rank %d: Data received from root process\n", rank);
//     }
// }

// // Step function
// void step() {
//     MPI_Request reqs[4];
//     MPI_Status stats[4];

//     printf("Rank %d: Starting ghost cell exchange\n", rank);

//     // Example ghost cell exchange (left-right neighbors)
//     if (rank > 0) {
//         MPI_Isend(h_local, ny_local, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &reqs[0]);
//         MPI_Irecv(h_local - ny_local, ny_local, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &reqs[1]);
//     }
//     if (rank < num_procs - 1) {
//         MPI_Isend(h_local + (nx_local - 1) * ny_local, ny_local, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &reqs[2]);
//         MPI_Irecv(h_local + nx_local * ny_local, ny_local, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &reqs[3]);
//     }

//     // Wait for communication to complete
//     MPI_Waitall(4, reqs, stats);
//     printf("Rank %d: Completed ghost cell exchange\n", rank);

//     // Update h_local, u_local, and v_local based on simulation logic
//     for (int i = 1; i < nx_local - 1; i++) {
//         for (int j = 1; j < ny_local - 1; j++) {
//             int idx = i * ny_local + j;

//             double dh_dx = (h_local[idx + ny_local] - h_local[idx - ny_local]) / 2.0;
//             double dh_dy = (h_local[idx + 1] - h_local[idx - 1]) / 2.0;
//             double du_dx = (u_local[idx + ny_local] - u_local[idx - ny_local]) / 2.0;
//             double dv_dy = (v_local[idx + 1] - v_local[idx - 1]) / 2.0;

//             h_local[idx] -= dt * (u_local[idx] * dh_dx + v_local[idx] * dh_dy);
//             u_local[idx] -= dt * (g * dh_dx + u_local[idx] * du_dx);
//             v_local[idx] -= dt * (g * dh_dy + v_local[idx] * dv_dy);
//         }
//     }
//     printf("Rank %d: Completed update step\n", rank);
// }

// // Transfer function
// void transfer(double *h_recv) {
//     printf("Rank %d: Starting data transfer to root\n", rank);
//     MPI_Gather(h_local, nx_local * ny_local, MPI_DOUBLE, h_recv, nx_local * ny_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     printf("Rank %d: Completed data transfer\n", rank);
// }

// // Free memory function
// void free_memory() {
//     printf("Rank %d: Freeing memory\n", rank);
//     free(h_local);
//     free(u_local);
//     free(v_local);
// }


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"

// Global variables
int nx, ny;
double *h, *u, *v, *dh, *du, *dv, *dh1, *du1, *dv1, *dh2, *du2, *dv2;
double H, g, dx, dy, dt;
int rank, num_procs;

// Initialize function
void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    h = h0;
    u = u0;
    v = v0;
    nx = nx_ / num_procs_;  // Divide nx by the number of processors
    ny = ny_;

    dh = (double *)calloc(nx * ny, sizeof(double));
    du = (double *)calloc(nx * ny, sizeof(double));
    dv = (double *)calloc(nx * ny, sizeof(double));

    dh1 = (double *)calloc(nx * ny, sizeof(double));
    du1 = (double *)calloc(nx * ny, sizeof(double));
    dv1 = (double *)calloc(nx * ny, sizeof(double));

    dh2 = (double *)calloc(nx * ny, sizeof(double));
    du2 = (double *)calloc(nx * ny, sizeof(double));
    dv2 = (double *)calloc(nx * ny, sizeof(double));

    H = H_;
    g = g_;
    dx = length_ / (nx * num_procs_);  // Account for global domain size
    dy = width_ / ny;
    dt = dt_;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
}

// Helper macro to access 2D arrays stored in 1D format
#define IDX(i, j) ((i) * ny + (j))

// Ghost cell exchange
void exchange_ghost_cells()
{
    MPI_Request requests[4] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};

    // Allocate boundary buffers to handle ghost cell data
    double *left_buffer = (double *)malloc(ny * sizeof(double));
    double *right_buffer = (double *)malloc(ny * sizeof(double));

    if (left_buffer == NULL || right_buffer == NULL) {
        fprintf(stderr, "Rank %d: Failed to allocate ghost cell buffers.\n", rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // Send right boundary, receive left boundary
    if (rank < num_procs - 1) {
        printf("Rank %d sending right boundary to %d\n", rank, rank + 1);
        MPI_Isend(&h[IDX(nx - 1, 0)], ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &requests[0]);
    }
    if (rank > 0) {
        printf("Rank %d receiving left boundary from %d\n", rank, rank - 1);
        MPI_Irecv(left_buffer, ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &requests[1]);
    }

    // Send left boundary, receive right boundary
    if (rank > 0) {
        printf("Rank %d sending left boundary to %d\n", rank, rank - 1);
        MPI_Isend(&h[IDX(0, 0)], ny, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &requests[2]);
    }
    if (rank < num_procs - 1) {
        printf("Rank %d receiving right boundary from %d\n", rank, rank + 1);
        MPI_Irecv(right_buffer, ny, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &requests[3]);
    }

    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

    // Check if received buffers are populated
    if (rank > 0) {
        for (int j = 0; j < ny; j++) {
            h[IDX(-1, j)] = left_buffer[j];
        }
    }
    if (rank < num_procs - 1) {
        for (int j = 0; j < ny; j++) {
            h[IDX(nx, j)] = right_buffer[j];
        }
    }

    // Free buffers
    free(left_buffer);
    free(right_buffer);
}



void compute_dh()
{
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            dh[IDX(i, j)] = -H * (u[IDX(i + 1, j)] - u[IDX(i - 1, j)]) / (2.0 * dx) - (v[IDX(i, j + 1)] - v[IDX(i, j - 1)]) / (2.0 * dy);
        }
    }
}

void compute_du()
{
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            du[IDX(i, j)] = -g * (h[IDX(i + 1, j)] - h[IDX(i - 1, j)]) / (2.0 * dx);
        }
    }
}

void compute_dv()
{
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            dv[IDX(i, j)] = -g * (h[IDX(i, j + 1)] - h[IDX(i, j - 1)]) / (2.0 * dy);
        }
    }
}

void multistep(double a1, double a2, double a3)
{
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            h[IDX(i, j)] += (a1 * dh[IDX(i, j)] + a2 * dh1[IDX(i, j)] + a3 * dh2[IDX(i, j)]) * dt;
            u[IDX(i, j)] += (a1 * du[IDX(i, j)] + a2 * du1[IDX(i, j)] + a3 * du2[IDX(i, j)]) * dt;
            v[IDX(i, j)] += (a1 * dv[IDX(i, j)] + a2 * dv1[IDX(i, j)] + a3 * dv2[IDX(i, j)]) * dt;
        }
    }
}

// Swap buffers for multistep method
void swap_buffers()
{
    double *tmp;

    tmp = dh2;
    dh2 = dh1;
    dh1 = dh;
    dh = tmp;

    tmp = du2;
    du2 = du1;
    du1 = du;
    du = tmp;

    tmp = dv2;
    dv2 = dv1;
    dv1 = dv;
    dv = tmp;
}

// Simulation step
void step()
{
    static int t = 0;
    exchange_ghost_cells();
    compute_dh();
    compute_du();
    compute_dv();

    double a1, a2 = 0.0, a3 = 0.0;  // Initialize a2, a3 to avoid undefined behavior

    if (t == 0) {
        a1 = 1.0;
    } else if (t == 1) {
        a1 = 3.0 / 2.0;
        a2 = -1.0 / 2.0;
    } else {
        a1 = 23.0 / 12.0;
        a2 = -16.0 / 12.0;
        a3 = 5.0 / 12.0;
    }

    multistep(a1, a2, a3);
    swap_buffers();
    t++;
}

// Data transfer to root process
void transfer(double *h_recv)
{
    int local_size = nx * ny;
    double *global_h = nullptr;

    if (rank == 0) {
        global_h = (double *)malloc(nx * ny * num_procs * sizeof(double));
    }

    MPI_Gather(h, local_size, MPI_DOUBLE, global_h, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Process global_h as needed
        free(global_h);
    }
}

// Free allocated memory
void free_memory()
{
    free(dh);
    free(du);
    free(dv);
    free(dh1);
    free(du1);
    free(dv1);
    free(dh2);
    free(du2);
    free(dv2);
}
