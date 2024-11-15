#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

#define ny local_ny
#include "../common/common.hpp"
#include "../common/solver.hpp"

int rank, num_procs;
int local_nx, local_ny;
int global_nx, global_ny;

double *h, *u, *v;
double *dh, *du, *dv;
double *dh1, *du1, *dv1;
double *dh2, *du2, *dv2;
double *send_buffer, *recv_buffer;

double H, g, dx, dy, dt;
int t = 0;

void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
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

    h = (double*)calloc((local_nx+2) * (local_ny+2), sizeof(double));
    u = (double*)calloc((local_nx+2) * (local_ny+2), sizeof(double));
    v = (double*)calloc((local_nx+2) * (local_ny+2 ), sizeof(double));
    
    dh = (double*)calloc((local_nx) * (local_ny), sizeof(double));
    du = (double*)calloc((local_nx) * (local_ny), sizeof(double));
    dv = (double*)calloc((local_nx) * (local_ny), sizeof(double));
    
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
        sendcounts[i] = nx_i * (global_ny + 1);
        displs[i] = offset;
        offset += sendcounts[i];
    }

    if (rank == 0) {
        MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(u0, sendcounts, displs, MPI_DOUBLE, u, local_nx * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, u, local_nx * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    delete[] sendcounts;
    delete[] displs;
}

void exchange_ghost_cells() {
    MPI_Status status;
    int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
    int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;

    // Exchange h ghost cells
    // Send rightmost column to right neighbor
    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = h(local_nx-1, j);  // Last real column
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 0,
                 recv_buffer, local_ny, MPI_DOUBLE, left, 0,
                 MPI_COMM_WORLD, &status);
    // Receive into leftmost ghost column
    for (int j = 0; j < local_ny; j++) {
        h(0, j) = recv_buffer[j];
    }

    // Exchange back (left to right)
    for (int j = 0; j < local_ny; j++) {
        send_buffer[j] = h(1, j);  // First real column
    }
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 1,
                 recv_buffer, local_ny, MPI_DOUBLE, right, 1,
                 MPI_COMM_WORLD, &status);
    // Receive into rightmost ghost column  
    for (int j = 0; j < local_ny; j++) {
        h(local_nx, j) = recv_buffer[j];
    }

}

void compute_derivatives()
{
    for (int i = 0; i < local_nx; i++)
    {
        for (int j = 0; j < local_ny+1; j++)
        {
            double dhdx = (h(i+1, j) - h(i, j)) / dx;
            double dhdy = (h(i, j + 1) - h(i, j)) / dy;

            dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
            du(i, j) = -g * dhdx;
            dv(i, j) = -g * dhdy;
        }
    }
}

void step()
{
    exchange_ghost_cells();
    compute_derivatives();

    double a1, a2, a3;
    if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
    else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
    else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


    for (int i = 0; i < local_nx; i++)
    {
        for (int j = 0; j < local_ny; j++)
        {
            h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
            u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
            v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
        }
    }

    double *tmp;
    tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
    tmp = du2; du2 = du1; du1 = du; du = tmp;
    tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

    t++;
}

void transfer(double *h_recv)
{
    int *recvcounts = new int[num_procs];
    int *displs = new int[num_procs];
    int offset = 0;

    for (int i = 0; i < num_procs; i++) {
        int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
        recvcounts[i] = nx_i * (global_ny + 1);
        displs[i] = offset;
        offset += recvcounts[i];
    }

    MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    delete[] recvcounts;
    delete[] displs;
}

void free_memory() {
    // Track what's been freed
    static bool freed = false;
    if (freed) return;
    
    if (h) { free(h); h = nullptr; }
    if (u) { free(u); u = nullptr; }
    if (v) { free(v); v = nullptr; }
    
    if (dh) { free(dh); dh = nullptr; }
    if (du) { free(du); du = nullptr; }
    if (dv) { free(dv); dv = nullptr; }
    
    if (dh1) { free(dh1); dh1 = nullptr; }
    if (du1) { free(du1); du1 = nullptr; }
    if (dv1) { free(dv1); dv1 = nullptr; }
    
    if (dh2) { free(dh2); dh2 = nullptr; }
    if (du2) { free(du2); du2 = nullptr; }
    if (dv2) { free(dv2); dv2 = nullptr; }
    
    if (send_buffer) { free(send_buffer); send_buffer = nullptr; }
    if (recv_buffer) { free(recv_buffer); recv_buffer = nullptr; }
    
    freed = true;
}