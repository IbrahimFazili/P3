// // max error 8.881784197001252e-15
// // acs378@nid001005:~/P3> python ./utils/timing.py gpu build/gpu
// // Timing gpu size 100
// // Total time: 1.2717420000000002
// // Timing gpu size 500
// // Total time: 1.308444
// // Timing gpu size 1000
// // Total time: 3.187059
// // Timing gpu size 2000
// // Total time: 9.214706999999999
// // Timing gpu size 3000
// // Total time: 19.502035
// // Timing gpu size 4000
// // Total time: 33.185929
// // Timing gpu size 5000
// // Total time: 53.000116
// // Timing gpu size 6000
// // Total time: 72.269965
// // Timing gpu size 7000
// // Total time: 98.306725
// // Timing gpu size 8000
// // Total time: 125.49425200000002
// // Timing gpu size 9000
// // Total time: 1.841289
// // Timing gpu size 10000
// // Total time: 192.431113
// #include <cuda.h>
// #include <cuda_runtime.h>
// #include <math.h>
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// // vars for grid size
// int nx, ny;
// double H, g, dx, dy, dt;

// // device ptrs for fields and derivs
// double *d_h, *d_u, *d_v, *d_dh, *d_du, *d_dv, *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
// int t = 0;

// // kernel to compute dh on GPU
// __global__ void compute_dh(double *d_dh, double *d_u, double *d_v, double H, double dx, double dy, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
    
//     if (i < nx && j < ny) {
//         double du_dx = (d_u[(i + 1) * ny + j] - d_u[i * ny + j]) / dx;
//         double dv_dy = (d_v[i * (ny + 1) + j + 1] - d_v[i * (ny + 1) + j]) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }
// }

// // kernel to compute du on GPU
// __global__ void compute_du(double *d_du, double *d_h, double g, double dx, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
    
//     if (i < nx && j < ny) {
//         double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - d_h[i * (ny + 1) + j]) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }
// }

// // kernel to compute dv on GPU
// __global__ void compute_dv(double *d_dv, double *d_h, double g, double dy, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
    
//     if (i < nx && j < ny) {
//         double dh_dy = (d_h[i * (ny + 1) + j + 1] - d_h[i * (ny + 1) + j]) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
// }

// // kernel to perform multistep update for h, u, v
// __global__ void multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                           double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                           double a1, double a2, double a3, double dt, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
    
//     if (i < nx && j < ny) {
//         d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
//         d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//         d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//     }
// }

// // init fn to allocate memory on GPU and copy initial data
// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     // Set grid parameters
//     nx = nx_;
//     ny = ny_;
//     H = H_;
//     g = g_;
//     dx = length_ / nx;
//     dy = width_ / ny;
//     dt = dt_;

//     // allocate device memory
//     cudaMalloc(&d_h, (nx + 1) * (ny + 1) * sizeof(double));
//     cudaMalloc(&d_u, nx * ny * sizeof(double));
//     cudaMalloc(&d_v, (nx + 1) * (ny + 1) * sizeof(double));

//     cudaMalloc(&d_dh, nx * ny * sizeof(double));
//     cudaMalloc(&d_du, nx * ny * sizeof(double));
//     cudaMalloc(&d_dv, nx * ny * sizeof(double));

//     cudaMalloc(&d_dh1, nx * ny * sizeof(double));
//     cudaMalloc(&d_du1, nx * ny * sizeof(double));
//     cudaMalloc(&d_dv1, nx * ny * sizeof(double));

//     cudaMalloc(&d_dh2, nx * ny * sizeof(double));
//     cudaMalloc(&d_du2, nx * ny * sizeof(double));
//     cudaMalloc(&d_dv2, nx * ny * sizeof(double));

//     // cp init data to GPU
//     cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_u, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_v, v0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
// }

// // step function to perform single time step on GPU
// void step()
// {
//     // set up thread block and grid dimensions
//     dim3 blockSize(16, 16);
//     dim3 gridSize((nx + blockSize.x - 1) / blockSize.x, (ny + blockSize.y - 1) / blockSize.y);

//     // compute coefficients for multistep method
//     double a1, a2, a3;
//     if (t == 0) {
//         a1 = 1.0;
//         a2 = 0.0;
//         a3 = 0.0;
//     } else if (t == 1) {
//         a1 = 3.0 / 2.0;
//         a2 = -1.0 / 2.0;
//         a3 = 0.0;
//     } else {
//         a1 = 23.0 / 12.0;
//         a2 = -16.0 / 12.0;
//         a3 = 5.0 / 12.0;
//     }

//     // launch kernels to compute derivatives
//     compute_dh<<<gridSize, blockSize>>>(d_dh, d_u, d_v, H, dx, dy, nx, ny);
//     compute_du<<<gridSize, blockSize>>>(d_du, d_h, g, dx, nx, ny);
//     compute_dv<<<gridSize, blockSize>>>(d_dv, d_h, g, dy, nx, ny);

//     // launch kernel to update h, u, v w/ multistep method
//     multistep<<<gridSize, blockSize>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2, a1, a2, a3, dt, nx, ny);

//     // swap derivative buffers
//     double *tmp;
//     tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
//     tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
//     tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;

//     t++;
// }

// // transfer fn to copy the h field back to host
// void transfer(double *h_host)
// {
//     cudaMemcpy(h_host, d_h, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyDeviceToHost);
// }

// // free GPU memory
// void free_memory()
// {
//     cudaFree(d_h);
//     cudaFree(d_u);
//     cudaFree(d_v);

//     cudaFree(d_dh);
//     cudaFree(d_du);
//     cudaFree(d_dv);

//     cudaFree(d_dh1);
//     cudaFree(d_du1);
//     cudaFree(d_dv1);

//     cudaFree(d_dh2);
//     cudaFree(d_du2);
//     cudaFree(d_dv2);
// }

#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <cstdio>     // fprintf
#include <cstdlib>    // exit
#include "../common/common.hpp"
#include "../common/solver.hpp"

// vars for grid size
int nx, ny;
double H, g, dx, dy, dt;

// dev ptrs for fields and derivs
double *d_h, *d_u, *d_v, *d_dh, *d_du, *d_dv, *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
int t = 0;

// macro to check for cuda errors
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error in %s at line %d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(err); \
        } \
    } while (0)

// kernel to compute dh on gpu w boundary checks
__global__ void compute_dh(double *d_dh, double *d_u, double *d_v, double H, double dx, double dy, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx - 1 && j < ny - 1) { // avoid out-of-bounds access
        double du_dx = (d_u[(i + 1) * ny + j] - d_u[i * ny + j]) / dx;
        double dv_dy = (d_v[i * (ny + 1) + j + 1] - d_v[i * (ny + 1) + j]) / dy;
        d_dh[i * ny + j] = -H * (du_dx + dv_dy);
    }
}

// kernel to compute du on gpu w boundary checks
__global__ void compute_du(double *d_du, double *d_h, double g, double dx, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx - 1 && j < ny) { // avoid out-of-bounds access
        double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - d_h[i * (ny + 1) + j]) / dx;
        d_du[i * ny + j] = -g * dh_dx;
    }
}

// kernel to compute dv on gpu w boundary checks
__global__ void compute_dv(double *d_dv, double *d_h, double g, double dy, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx && j < ny - 1) { // avoid out-of-bounds access
        double dh_dy = (d_h[i * (ny + 1) + j + 1] - d_h[i * (ny + 1) + j]) / dy;
        d_dv[i * ny + j] = -g * dh_dy;
    }
}

// kernel to perform multistep update for h, u, v w additional boundary checks
__global__ void multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
                          double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
                          double a1, double a2, double a3, double dt, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx && j < ny) { // Ensure within bounds
        d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
        
        // Add boundary checks for d_u and d_v accesses
        if (i + 1 < nx) { // Avoid out-of-bounds access for d_u
            d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
        }
        if (j + 1 < ny) { // Avoid out-of-bounds access for d_v
            d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
        }
    }
}

// init fn to allocate mem on gpu and cp initial data 
void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    // set grid parameters
    nx = nx_;
    ny = ny_;
    H = H_;
    g = g_;
    dx = length_ / nx;
    dy = width_ / ny;
    dt = dt_;

    // alloc device mem
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

    // cp init data to gpu
    cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_u, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, v0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
}

// step fn to perform single time step on gpu 
void step()
{
    // set up thread block and grid dimensions
    dim3 blockSize(16, 16);
    dim3 gridSize((nx + blockSize.x - 1) / blockSize.x, (ny + blockSize.y - 1) / blockSize.y);

    // compute coefficients for multistep method
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

    // launch kernels to compute derivs
    compute_dh<<<gridSize, blockSize>>>(d_dh, d_u, d_v, H, dx, dy, nx, ny);
    cudaDeviceSynchronize();
    
    compute_du<<<gridSize, blockSize>>>(d_du, d_h, g, dx, nx, ny);
    cudaDeviceSynchronize();

    compute_dv<<<gridSize, blockSize>>>(d_dv, d_h, g, dy, nx, ny);
    cudaDeviceSynchronize();

    // launch kernel to update h, u, v with multistep method
    multistep<<<gridSize, blockSize>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2, a1, a2, a3, dt, nx, ny);
    cudaDeviceSynchronize();

    // swap deriv buffers
    double *tmp;
    tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
    tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
    tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;

    t++;
}

// transfer fn to copy h field back to host 
void transfer(double *h_host)
{
    cudaMemcpy(h_host, d_h, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyDeviceToHost);
}

// free gpu mem 
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
