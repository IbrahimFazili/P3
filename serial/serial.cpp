#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"

// vars for grid size
int nx, ny;

// arrays to hold field vals and derivs
double *h, *u, *v, *dh, *du, *dv, *dh1, *du1, *dv1, *dh2, *du2, *dv2;
double H, g, dx, dy, dt;

int t = 0;

/**
 * init fn
 */
void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    // ptrs to the arrays passed in
    h = h0;
    u = u0;
    v = v0;

    nx = nx_;
    ny = ny_;

    // alloc memory for derivates
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

    dx = length_ / nx;
    dy = width_ / ny;
    dt = dt_;
}

/**
 * compute derivs of the fields in a single pass (for cache efficiency)
 */
void compute_derivatives()
{
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            double dhdx = (h(i + 1, j) - h(i, j)) / dx;
            double dhdy = (h(i, j + 1) - h(i, j)) / dy;

            dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
            du(i, j) = -g * dhdx;
            dv(i, j) = -g * dhdy;
        }
    }
}

/**
 * compute the next time step using multistep method
 */
void multistep(double a1, double a2, double a3)
{
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
            u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
            v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
        }
    }
}

/**
 * step fm
 */
void step()
{
    // compute ghost cells
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) h(nx, j) = h(0, j);
    for (int i = 0; i < nx; i++) h(i, ny) = h(i, 0);

    // computs derivs
    compute_derivatives();

    // set coefficients for multistep method based on time step
    double a1, a2, a3;
    if (t == 0) {
        a1 = 1.0; a2 = 0.0; a3 = 0.0;
    } else if (t == 1) {
        a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0;
    } else {
        a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0;
    }

    // compute next time step
    multistep(a1, a2, a3);

    // update boundaries
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) u(0, j) = u(nx, j);
    for (int i = 0; i < nx; i++) v(i, 0) = v(i, ny);

    // swap buffers for multistep method
    double *tmp;
    tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
    tmp = du2; du2 = du1; du1 = du; du = tmp;
    tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

    t++;
}

/**
 * This is your transfer function! Since everything is running on the same node,
 * you don't need to do anything here.
 */
void transfer(double *h)
{
    return;
}

/**
 * This is your finalization function! Free whatever memory you've allocated.
 */
void free_memory()
{
    free(dh); free(du); free(dv);
    free(dh1); free(du1); free(dv1);
    free(dh2); free(du2); free(dv2);
}