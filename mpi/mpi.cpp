#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"
#define idx(i, j) ((i) * (local_ny) + (j))

int rank, num_procs;
int local_nx, local_ny;
int global_nx, global_ny;
int start_x, end_x;

double *h, *u, *v;
double *dh, *du, *dv;
double *dh1, *du1, *dv1;
double *dh2, *du2, *dv2;
double *send_buffer, *recv_buffer;

double H, g, dx, dy, dt;
int t = 0;

/**
 * This is your initialization function! It is very similar to the one in
 * serial.cpp, but with some difference. Firstly, only the process with rank 0
 * is going to actually generate the initial conditions h0, u0, and v0, so all
 * other processes are going to get nullptrs. Therefore, you'll need to find some
 * way to scatter the initial conditions to all processes. Secondly, now the
 * rank and num_procs arguments are passed to the function, so you can use them
 * to determine which rank the node running this process has, and how many
 * processes are running in total. This is useful to determine which part of the
 * domain each process is going to be responsible for.
 */
void init(double *h0, double *u0, double *v0, double length_, double width_, 
          int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    rank = rank_;
    num_procs = num_procs_;
    global_nx = nx_;
    global_ny = ny_;
    
    local_nx = global_nx / num_procs;
    if (rank < global_nx % num_procs) {
        local_nx++;
    }
    local_ny = global_ny;
    
    start_x = rank * (global_nx / num_procs);
    if (rank < global_nx % num_procs) {
        start_x += rank;
    } else {
        start_x += global_nx % num_procs;
    }
    end_x = start_x + local_nx;

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

    send_buffer = (double*)calloc(local_ny, sizeof(double));
    recv_buffer = (double*)calloc(local_ny, sizeof(double));

    H = H_;
    g = g_;
    dx = length_ / global_nx;
    dy = width_ / global_ny;
    dt = dt_;

    int *sendcounts = new int[num_procs];
    int *displs = new int[num_procs];
    int offset = 0;

    for (int i = 0; i < num_procs; i++) {
        int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
        sendcounts[i] = nx_i * global_ny;
        displs[i] = offset;
        offset += sendcounts[i];
    }

    if (rank == 0) {
        MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h + (local_ny + 2), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(u0, sendcounts, displs, MPI_DOUBLE, u + (local_ny + 2), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v + (local_ny + 2), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h + (local_ny + 2), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, u + (local_ny + 2), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v + (local_ny + 2), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    delete[] sendcounts;
    delete[] displs;

    double init_time = MPI_Wtime();
    init_time = MPI_Wtime() - init_time;
    printf("Initialization time for rank %d: %.6f\n", rank, init_time);
}

void exchange_ghost_cells()
{
    MPI_Status status;
    int left = (rank == 0) ? num_procs - 1 : rank - 1;
    int right = (rank == num_procs - 1) ? 0 : rank + 1;

    // Exchange h ghost cells
    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = h[idx(local_nx, j+1)];
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 0,
                 recv_buffer, local_ny, MPI_DOUBLE, left, 0,
                 MPI_COMM_WORLD, &status);
    for (int j = 0; j < local_ny; j++) {
        h[idx(0, j+1)] = recv_buffer[j];
    }

    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = h[idx(1, j+1)];
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 1,
                 recv_buffer, local_ny, MPI_DOUBLE, right, 1,
                 MPI_COMM_WORLD, &status);
    for (int j = 0; j < local_ny; j++) {
        h[idx(local_nx+1, j+1)] = recv_buffer[j];
    }

    // Exchange u ghost cells
    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = u[idx(local_nx, j+1)];
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 2,
                 recv_buffer, local_ny, MPI_DOUBLE, left, 2,
                 MPI_COMM_WORLD, &status);
    for (int j = 0; j < local_ny; j++) {
        u[idx(0, j+1)] = recv_buffer[j];
    }

    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = u[idx(1, j+1)];
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 3,
                 recv_buffer, local_ny, MPI_DOUBLE, right, 3,
                 MPI_COMM_WORLD, &status);
    for (int j = 0; j < local_ny; j++) {
        u[idx(local_nx+1, j+1)] = recv_buffer[j];
    }

    // Exchange v ghost cells
    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = v[idx(local_nx, j+1)];
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 4,
                 recv_buffer, local_ny, MPI_DOUBLE, left, 4,
                 MPI_COMM_WORLD, &status);
    for (int j = 0; j < local_ny; j++) {
        v[idx(0, j+1)] = recv_buffer[j];
    }

    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = v[idx(1, j+1)];
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 5,
                 recv_buffer, local_ny, MPI_DOUBLE, right, 5,
                 MPI_COMM_WORLD, &status);
    for (int j = 0; j < local_ny; j++) {
        v[idx(local_nx+1, j+1)] = recv_buffer[j];
    }
}


void compute_derivatives()
{
    for (int i = 1; i <= local_nx; i++) {
        for (int j = 1; j <= local_ny; j++) {
            double dx_h = (h[idx(i+1,j)] - h[idx(i-1,j)]) / (2.0 * dx);
            double dy_h = (h[idx(i,j+1)] - h[idx(i,j-1)]) / (2.0 * dy);
            
            double dx_u = (u[idx(i+1,j)] - u[idx(i-1,j)]) / (2.0 * dx);
            double dy_v = (v[idx(i,j+1)] - v[idx(i,j-1)]) / (2.0 * dy);

            dh[idx(i-1,j-1)] = -H * (dx_u + dy_v);
            du[idx(i-1,j-1)] = -g * dx_h;
            dv[idx(i-1,j-1)] = -g * dy_h;
        }
    }
}

/**
 * This is your step function! It is very similar to the one in serial.cpp, but
 * now the domain is divided among the processes, so you'll need to find some
 * way to communicate the ghost cells between processes.
 */
void step()
{
    exchange_ghost_cells();
    compute_derivatives();

    double a1, a2, a3;
    if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
    else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
    else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }

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

    double *tmp;
    tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
    tmp = du2; du2 = du1; du1 = du; du = tmp;
    tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

    t++;
}

/**
 * This is your transfer function! Similar to what you did in gpu.cu, you'll
 * need to get the data from the computers you're working on (there it was
 * the GPU, now its a bunch of CPU nodes), and send them all back to the process
 * which is actually running the main function (then it was the CPU, not it's
 * the node with rank 0).
 */
void transfer(double *h_recv)
{
    int *recvcounts = new int[num_procs];
    int *displs = new int[num_procs];
    int offset = 0;

    for (int i = 0; i < num_procs; i++) {
        int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
        recvcounts[i] = nx_i * local_ny;
        displs[i] = offset;
        offset += recvcounts[i];
    }

    MPI_Gatherv(h + (local_ny + 2), local_nx * local_ny, MPI_DOUBLE, 
                h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    delete[] recvcounts;
    delete[] displs;
}


/**
 * This is your finalization function! Since different nodes are going to be
 * initializing different chunks of memory, make sure to check which node
 * is running the code before you free some memory you haven't allocated, or
 * that you've actually freed memory that you have.
 */
void free_memory()
{
    free(h); free(u); free(v);
    free(dh); free(du); free(dv);
    free(dh1); free(du1); free(dv1);
    free(dh2); free(du2); free(dv2);
    free(send_buffer);
    free(recv_buffer);
}