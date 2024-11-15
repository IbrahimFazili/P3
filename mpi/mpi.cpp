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

    h = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
    u = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
    v = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
    dh = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
    du = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
    dv = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
    dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
    dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
    du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
    dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

    send_buffer = (double*)calloc(local_ny+2, sizeof(double));
    recv_buffer = (double*)calloc(local_ny+2, sizeof(double));

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

void exchange_ghost_cells()
{
    MPI_Status status;
    int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
    int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;

    std::memcpy(send_buffer, h + local_nx * local_ny, local_ny * sizeof(double));

    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 0, recv_buffer, local_ny, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
    
    for (int j = 1; j <= local_ny; j++) {
        h(local_nx + 1, j) = recv_buffer[j - 1];
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

void free_memory()
{
    free(h); free(u); free(v);
    free(dh); free(du); free(dv);
    free(dh1); free(du1); free(dv1);
    free(dh2); free(du2); free(dv2);
    free(send_buffer);
    free(recv_buffer);
}