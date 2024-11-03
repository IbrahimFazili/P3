#include <cuda.h>
#include <cuda_runtime.h>

#include <math.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"

/**
 * This is your initialization function! We pass in h0, u0, and v0, which are
 * your initial height, u velocity, and v velocity fields. You should send these
 * grids to the GPU so you can do work on them there, and also these other fields.
 * Here, length and width are the length and width of the domain, and nx and ny are
 * the number of grid points in the x and y directions. H is the height of the water
 * column, g is the acceleration due to gravity, and dt is the time step size.
 * The rank and num_procs variables are unused here, but you will need them
 * when doing the MPI version.
 */

// Here we hold the number of cells we have in the x and y directions
int nx, ny;

double *h, *u, *v, *d_dh, *d_du, *d_dv;
double *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
double H, g, dx, dy, dt;
double *cpu_h;

void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    // Set the grid size
    nx = nx_;
    ny = ny_;

    // Allocate memory on the device
    cudaMalloc((void**)&h, nx * ny * sizeof(double));
    cudaMalloc((void**)&u, nx * ny * sizeof(double));
    cudaMalloc((void**)&v, nx * ny * sizeof(double));

    cudaMalloc((void**)&d_dh, nx * ny * sizeof(double));
    cudaMalloc((void**)&d_du, nx * ny * sizeof(double));
    cudaMalloc((void**)&d_dv, nx * ny * sizeof(double));

    cudaMalloc((void**)&d_dh1, nx * ny * sizeof(double));
    cudaMalloc((void**)&d_du1, nx * ny * sizeof(double));
    cudaMalloc((void**)&d_dv1, nx * ny * sizeof(double));

    cudaMalloc((void**)&d_dh2, nx * ny * sizeof(double));
    cudaMalloc((void**)&d_du2, nx * ny * sizeof(double));
    cudaMalloc((void**)&d_dv2, nx * ny * sizeof(double));

    // Copy initial conditions to the device @TODO: change

    cudaMemcpy(d_dh, h0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_du, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dv, v0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(d_dh1, h0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_du1, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dv1, v0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(d_dh2, h0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_du2, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dv2, v0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);

    cpu_h = (double*)calloc(nx * ny, sizeof(double));

    H = H_;
    g = g_;

    dx = length_ / nx;
    dy = width_ / ny;

    dt = dt_;
}

/**
 * This is your step function! Here, you will actually numerically solve the shallow
 * water equations. You should update the h, u, and v fields to be the solution after
 * one time step has passed.
 */
void step()
{
    // @TODO: Your code here
}

__global__ void compute_ghost_horizontal(double *h, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    // @TOOD strides? for loop for everything below

    // we don't want to go out of bounds
    if (i < ny) {
        h(nx, i) = h(0, i);
    }
}

__global__ void compute_ghost_vertical(double *h, int nx, int ny) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (j < nx) {
        h(j, ny) = h(j, 0);
    }
}

__global__ void compute_boundaries_horizontal(double *u, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < ny) {
        u(0, i) = u(nx, i);
    }
}

__global__ void compute_boundaries_vertical(double *v, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < nx) {
        v(i, 0) = v(i, ny);
    }
}

__global__ void compute_dh(
    double *dh, double *du, double *dv, 
    double *u, double *v, double H, double g, 
    int nx, int ny)  {
    // will require loop
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int stride = blockDim.x;

    for (; i < nx; i += stride) {
        for (; j < ny; j += stride) {
            // @TODO: Your code here
            dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
            du(i, j) = -g * dh_dx(i, j);
            dv(i, j) = -g * dh_dy(i, j);
        }
    }
}

/**
 * This is your transfer function! You should copy the h field back to the host
 * so that the CPU can check the results of your computation.
 */
void transfer(double *h_host)
{
    // @TODO: Your code here
}

/**
 * This is your finalization function! You should free all of the memory that you
 * allocated on the GPU here.
 */
void free_memory()
{
    cudaFree(h);
    cudaFree(u);
    cudaFree(v);

    cudaFree(d_dh);
    cudaFree(d_du);
    cudaFree(d_dv);

    cudaFree(d_dh1);
    cudaFree(d_du1);
    cudaFree(d_dv1);
    
    cudaFree(d_dh2);
    cudaFree(d_du2);
    cudaFree(d_dv2);

    free(cpu_h);
}