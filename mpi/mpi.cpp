#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"
// Add this near the top of your file, after the includes
#define idx(i,j) ((i) * (local_ny + 2) + (j))

// Global variables for domain decomposition
int rank, num_procs;
int local_nx, local_ny;
int global_nx, global_ny;
int start_x, end_x;

// Fields and derivatives
double *h, *u, *v;
double *dh, *du, *dv;
double *dh1, *du1, *dv1;
double *dh2, *du2, *dv2;
double *send_buffer, *recv_buffer;

// Physical parameters
double H, g, dx, dy, dt;
int t = 0;

void init(double *h0, double *u0, double *v0, double length_, double width_, 
         int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    rank = rank_;
    num_procs = num_procs_;
    global_nx = nx_;
    global_ny = ny_;
    
    // Calculate local domain size
    local_nx = global_nx / num_procs;
    if (rank < global_nx % num_procs) {
        local_nx++;
    }
    local_ny = global_ny;
    
    // Calculate starting and ending indices
    start_x = rank * (global_nx / num_procs);
    if (rank < global_nx % num_procs) {
        start_x += rank;
    } else {
        start_x += global_nx % num_procs;
    }
    end_x = start_x + local_nx;

    // Allocate local arrays with ghost cells
    h = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
    u = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
    v = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
    
    dh = (double*)calloc(local_nx * local_ny, sizeof(double));
    du = (double*)calloc(local_nx * local_ny, sizeof(double));
    dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
    dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
    dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
    du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
    dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

    // Buffers for MPI communication
    send_buffer = (double*)calloc(local_ny, sizeof(double));
    recv_buffer = (double*)calloc(local_ny, sizeof(double));

    // Initialize physical parameters
    H = H_;
    g = g_;
    dx = length_ / global_nx;
    dy = width_ / global_ny;
    dt = dt_;

    // Scatter initial conditions from rank 0
    if (rank == 0) {
        // Copy initial data to local arrays
        for (int i = 0; i < local_nx; i++) {
            for (int j = 0; j < local_ny; j++) {
                h[idx(i+1,j+1)] = h0[idx(i,j)];
                u[idx(i+1,j+1)] = u0[idx(i,j)];
                v[idx(i+1,j+1)] = v0[idx(i,j)];
            }
        }
    }

    // Broadcast initial conditions to all processes
    MPI_Bcast(h, (local_nx + 2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(u, (local_nx + 2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(v, (local_nx + 2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void exchange_ghost_cells()
{
    MPI_Status status;
    int left = (rank == 0) ? num_procs - 1 : rank - 1;
    int right = (rank == num_procs - 1) ? 0 : rank + 1;

    // Send right boundary, receive left ghost cells
    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = h[idx(local_nx, j+1)];
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 0,
                 recv_buffer, local_ny, MPI_DOUBLE, left, 0,
                 MPI_COMM_WORLD, &status);
    for (int j = 0; j < local_ny; j++) {
        h[idx(0, j+1)] = recv_buffer[j];
    }

    // Send left boundary, receive right ghost cells
    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = h[idx(1, j+1)];
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 1,
                 recv_buffer, local_ny, MPI_DOUBLE, right, 1,
                 MPI_COMM_WORLD, &status);
    for (int j = 0; j < local_ny; j++) {
        h[idx(local_nx+1, j+1)] = recv_buffer[j];
    }
}

void compute_derivatives()
{
    for (int i = 1; i <= local_nx; i++) {
        for (int j = 1; j <= local_ny; j++) {
            // Compute height derivatives
            double dx_h = (h[idx(i+1,j)] - h[idx(i-1,j)]) / (2.0 * dx);
            double dy_h = (h[idx(i,j+1)] - h[idx(i,j-1)]) / (2.0 * dy);
            
            // Compute velocity derivatives
            double dx_u = (u[idx(i+1,j)] - u[idx(i-1,j)]) / (2.0 * dx);
            double dy_v = (v[idx(i,j+1)] - v[idx(i,j-1)]) / (2.0 * dy);

            // Store derivatives
            dh[idx(i-1,j-1)] = -H * (dx_u + dy_v);
            du[idx(i-1,j-1)] = -g * dx_h;
            dv[idx(i-1,j-1)] = -g * dy_h;
        }
    }
}

void step()
{
    exchange_ghost_cells();
    compute_derivatives();

    // Set multistep coefficients
    double a1, a2, a3;
    if (t == 0) {
        a1 = 1.0;
        a2 = a3 = 0.0;
    } else if (t == 1) {
        a1 = 3.0 / 2.0;
        a2 = -1.0 / 2.0;
        a3 = 0.0;
    } else {
        a1 = 23.0 / 12.0;
        a2 = -16.0 / 12.0;
        a3 = 5.0 / 12.0;
    }

    // Update fields
    for (int i = 1; i <= local_nx; i++) {
        for (int j = 1; j <= local_ny; j++) {
            h[idx(i,j)] += (a1 * dh[idx(i-1,j-1)] + a2 * dh1[idx(i-1,j-1)] + 
                           a3 * dh2[idx(i-1,j-1)]) * dt;
            u[idx(i,j)] += (a1 * du[idx(i-1,j-1)] + a2 * du1[idx(i-1,j-1)] + 
                           a3 * du2[idx(i-1,j-1)]) * dt;
            v[idx(i,j)] += (a1 * dv[idx(i-1,j-1)] + a2 * dv1[idx(i-1,j-1)] + 
                           a3 * dv2[idx(i-1,j-1)]) * dt;
        }
    }

    // Swap derivative buffers
    double *tmp;
    tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
    tmp = du2; du2 = du1; du1 = du; du = tmp;
    tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

    t++;
}

void transfer(double *h_recv)
{
    if (rank == 0) {
        // Copy local data
        for (int i = 0; i < local_nx; i++) {
            for (int j = 0; j < local_ny; j++) {
                h_recv[i * local_ny + j] = h[idx(i+1,j+1)];
            }
        }

        // Receive data from other processes
        for (int p = 1; p < num_procs; p++) {
            int recv_size = (global_nx / num_procs) * local_ny;
            if (p < global_nx % num_procs) recv_size += local_ny;
            
            MPI_Recv(&h_recv[p * (global_nx / num_procs) * local_ny],
                    recv_size, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        // Send local data to rank 0
        MPI_Send(&h[local_ny + 1], local_nx * local_ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

void free_memory()
{
    free(h); free(u); free(v);
    free(dh); free(du); free(dv);
    free(dh1); free(du1); free(dv1);
    free(dh2); free(du2); free(dv2);
    free(send_buffer);
    free(recv_buffer);
}
