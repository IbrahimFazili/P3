// skeleton
// #include <cuda.h>
// #include <cuda_runtime.h>

// #include <math.h>

// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// /**
//  * This is your initialization function! We pass in h0, u0, and v0, which are
//  * your initial height, u velocity, and v velocity fields. You should send these
//  * grids to the GPU so you can do work on them there, and also these other fields.
//  * Here, length and width are the length and width of the domain, and nx and ny are
//  * the number of grid points in the x and y directions. H is the height of the water
//  * column, g is the acceleration due to gravity, and dt is the time step size.
//  * The rank and num_procs variables are unused here, but you will need them
//  * when doing the MPI version.
//  */
// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     // @TODO: your code here
// }

// /**
//  * This is your step function! Here, you will actually numerically solve the shallow
//  * water equations. You should update the h, u, and v fields to be the solution after
//  * one time step has passed.
//  */
// void step()
// {
//     // @TODO: Your code here
// }

// /**
//  * This is your transfer function! You should copy the h field back to the host
//  * so that the CPU can check the results of your computation.
//  */
// void transfer(double *h_host)
// {
//     // @TODO: Your code here
// }

// /**
//  * This is your finalization function! You should free all of the memory that you
//  * allocated on the GPU here.
//  */
// void free_memory()
// {
//     // @TODO: Your code here
// }


// ibrahim work
// #include <cuda.h>
// #include <cuda_runtime.h>

// #include <math.h>

// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// /**
//  * This is your initialization function! We pass in h0, u0, and v0, which are
//  * your initial height, u velocity, and v velocity fields. You should send these
//  * grids to the GPU so you can do work on them there, and also these other fields.
//  * Here, length and width are the length and width of the domain, and nx and ny are
//  * the number of grid points in the x and y directions. H is the height of the water
//  * column, g is the acceleration due to gravity, and dt is the time step size.
//  * The rank and num_procs variables are unused here, but you will need them
//  * when doing the MPI version.
//  */

// // Here we hold the number of cells we have in the x and y directions
// int nx, ny;

// double *h, *u, *v, *d_dh, *d_du, *d_dv;
// double *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
// double H, g, dx, dy, dt;
// double *cpu_h;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     // Set the grid size
//     nx = nx_;
//     ny = ny_;

//     // Allocate memory on the device
//     cudaMalloc((void**)&h, nx * ny * sizeof(double));
//     cudaMalloc((void**)&u, nx * ny * sizeof(double));
//     cudaMalloc((void**)&v, nx * ny * sizeof(double));

//     cudaMalloc((void**)&d_dh, nx * ny * sizeof(double));
//     cudaMalloc((void**)&d_du, nx * ny * sizeof(double));
//     cudaMalloc((void**)&d_dv, nx * ny * sizeof(double));

//     cudaMalloc((void**)&d_dh1, nx * ny * sizeof(double));
//     cudaMalloc((void**)&d_du1, nx * ny * sizeof(double));
//     cudaMalloc((void**)&d_dv1, nx * ny * sizeof(double));

//     cudaMalloc((void**)&d_dh2, nx * ny * sizeof(double));
//     cudaMalloc((void**)&d_du2, nx * ny * sizeof(double));
//     cudaMalloc((void**)&d_dv2, nx * ny * sizeof(double));

//     // Copy initial conditions to the device @TODO: change

//     cudaMemcpy(d_dh, h0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_du, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_dv, v0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);

//     cudaMemcpy(d_dh1, h0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_du1, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_dv1, v0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);

//     cudaMemcpy(d_dh2, h0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_du2, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_dv2, v0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);

//     cpu_h = (double*)calloc(nx * ny, sizeof(double));

//     H = H_;
//     g = g_;

//     dx = length_ / nx;
//     dy = width_ / ny;

//     dt = dt_;
// }

// /**
//  * This is your step function! Here, you will actually numerically solve the shallow
//  * water equations. You should update the h, u, and v fields to be the solution after
//  * one time step has passed.
//  */
// void step()
// {
//     // @TODO: Your code here
// }

// __global__ void compute_ghost_horizontal(double *h, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     // @TOOD strides? for loop for everything below

//     // we don't want to go out of bounds
//     if (i < ny) {
//         h(nx, i) = h(0, i);
//     }
// }

// __global__ void compute_ghost_vertical(double *h, int nx, int ny) {
//     int j = blockIdx.x * blockDim.x + threadIdx.x;

//     if (j < nx) {
//         h(j, ny) = h(j, 0);
//     }
// }

// __global__ void compute_boundaries_horizontal(double *u, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;

//     if (i < ny) {
//         u(0, i) = u(nx, i);
//     }
// }

// __global__ void compute_boundaries_vertical(double *v, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;

//     if (i < nx) {
//         v(i, 0) = v(i, ny);
//     }
// }

// __global__ void compute_dh(
//     double *dh, double *du, double *dv, 
//     double *u, double *v, double H, double g, 
//     int nx, int ny)  {
//     // will require loop
    
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
//     int stride = blockDim.x;

//     for (; i < nx; i += stride) {
//         for (; j < ny; j += stride) {
//             // @TODO: Your code here
//             dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
//             du(i, j) = -g * dh_dx(i, j);
//             dv(i, j) = -g * dh_dy(i, j);
//         }
//     }
// }

// /**
//  * This is your transfer function! You should copy the h field back to the host
//  * so that the CPU can check the results of your computation.
//  */
// void transfer(double *h_host)
// {
//     // @TODO: Your code here
// }

// /**
//  * This is your finalization function! You should free all of the memory that you
//  * allocated on the GPU here.
//  */
// void free_memory()
// {
//     cudaFree(h);
//     cudaFree(u);
//     cudaFree(v);

//     cudaFree(d_dh);
//     cudaFree(d_du);
//     cudaFree(d_dv);

//     cudaFree(d_dh1);
//     cudaFree(d_du1);
//     cudaFree(d_dv1);
    
//     cudaFree(d_dh2);
//     cudaFree(d_du2);
//     cudaFree(d_dv2);

//     free(cpu_h);
// }


// max error 8.881784197001252e-15
// acs378@nid001005:~/P3> python ./utils/timing.py gpu build/gpu
// Timing gpu size 100
// Total time: 1.2717420000000002
// Timing gpu size 500
// Total time: 1.308444
// Timing gpu size 1000
// Total time: 3.187059
// Timing gpu size 2000
// Total time: 9.214706999999999
// Timing gpu size 3000
// Total time: 19.502035
// Timing gpu size 4000
// Total time: 33.185929
// Timing gpu size 5000
// Total time: 53.000116
// Timing gpu size 6000
// Total time: 72.269965
// Timing gpu size 7000
// Total time: 98.306725
// Timing gpu size 8000
// Total time: 125.49425200000002
// Timing gpu size 9000
// Total time: 1.841289
// Timing gpu size 10000
// Total time: 192.431113
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include "../common/common.hpp"
#include "../common/solver.hpp"

// Variables for the grid size
int nx, ny;
double H, g, dx, dy, dt;

// Device pointers for fields and derivatives
double *d_h, *d_u, *d_v, *d_dh, *d_du, *d_dv, *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
int t = 0;

// Kernel to compute dh on the GPU
__global__ void compute_dh(double *d_dh, double *d_u, double *d_v, double H, double dx, double dy, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx && j < ny) {
        double du_dx = (d_u[(i + 1) * ny + j] - d_u[i * ny + j]) / dx;
        double dv_dy = (d_v[i * (ny + 1) + j + 1] - d_v[i * (ny + 1) + j]) / dy;
        d_dh[i * ny + j] = -H * (du_dx + dv_dy);
    }
}

// Kernel to compute du on the GPU
__global__ void compute_du(double *d_du, double *d_h, double g, double dx, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx && j < ny) {
        double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - d_h[i * (ny + 1) + j]) / dx;
        d_du[i * ny + j] = -g * dh_dx;
    }
}

// Kernel to compute dv on the GPU
__global__ void compute_dv(double *d_dv, double *d_h, double g, double dy, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx && j < ny) {
        double dh_dy = (d_h[i * (ny + 1) + j + 1] - d_h[i * (ny + 1) + j]) / dy;
        d_dv[i * ny + j] = -g * dh_dy;
    }
}

// Kernel to perform multistep update for h, u, and v
__global__ void multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
                          double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
                          double a1, double a2, double a3, double dt, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx && j < ny) {
        d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
        d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
        d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
    }
}

// Initialization function to allocate memory on the GPU and copy initial data
void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    // Set grid parameters
    nx = nx_;
    ny = ny_;
    H = H_;
    g = g_;
    dx = length_ / nx;
    dy = width_ / ny;
    dt = dt_;

    // Allocate device memory
    cudaMalloc(&d_h, (nx + 1) * (ny + 1) * sizeof(double));
    cudaMalloc(&d_u, nx * ny * sizeof(double));
    cudaMalloc(&d_v, (nx + 1) * (ny + 1) * sizeof(double));

    cudaMalloc(&d_dh, nx * ny * sizeof(double));
    cudaMalloc(&d_du, nx * ny * sizeof(double));
    cudaMalloc(&d_dv, nx * ny * sizeof(double));

    cudaMalloc(&d_dh1, nx * ny * sizeof(double));
    cudaMalloc(&d_du1, nx * ny * sizeof(double));
    cudaMalloc(&d_dv1, nx * ny * sizeof(double));

    cudaMalloc(&d_dh2, nx * ny * sizeof(double));
    cudaMalloc(&d_du2, nx * ny * sizeof(double));
    cudaMalloc(&d_dv2, nx * ny * sizeof(double));

    // Copy initial data to the GPU
    cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_u, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, v0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
}

// Step function to perform a single time step on the GPU
void step()
{
    // Set up thread block and grid dimensions
    dim3 blockSize(16, 16);
    dim3 gridSize((nx + blockSize.x - 1) / blockSize.x, (ny + blockSize.y - 1) / blockSize.y);

    // Compute the coefficients for the multistep method
    double a1, a2, a3;
    if (t == 0) {
        a1 = 1.0;
        a2 = 0.0;
        a3 = 0.0;
    } else if (t == 1) {
        a1 = 3.0 / 2.0;
        a2 = -1.0 / 2.0;
        a3 = 0.0;
    } else {
        a1 = 23.0 / 12.0;
        a2 = -16.0 / 12.0;
        a3 = 5.0 / 12.0;
    }

    // Launch kernels to compute derivatives
    compute_dh<<<gridSize, blockSize>>>(d_dh, d_u, d_v, H, dx, dy, nx, ny);
    compute_du<<<gridSize, blockSize>>>(d_du, d_h, g, dx, nx, ny);
    compute_dv<<<gridSize, blockSize>>>(d_dv, d_h, g, dy, nx, ny);

    // Launch kernel to update h, u, and v with the multistep method
    multistep<<<gridSize, blockSize>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2, a1, a2, a3, dt, nx, ny);

    // Swap derivative buffers
    double *tmp;
    tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
    tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
    tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;

    t++;
}

// Transfer function to copy the h field back to the host
void transfer(double *h_host)
{
    cudaMemcpy(h_host, d_h, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyDeviceToHost);
}

// Finalization function to free GPU memory
void free_memory()
{
    cudaFree(d_h);
    cudaFree(d_u);
    cudaFree(d_v);

    cudaFree(d_dh);
    cudaFree(d_du);
    cudaFree(d_dv);

    cudaFree(d_dh1);
    cudaFree(d_du1);
    cudaFree(d_dv1);

    cudaFree(d_dh2);
    cudaFree(d_du2);
    cudaFree(d_dv2);
}
