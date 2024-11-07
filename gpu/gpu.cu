#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <cstdio>     // fprintf
#include <cstdlib>    // exit
#include "../common/common.hpp"
#include "../common/solver.hpp"
#include <cublas_v2.h>

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
template <const uint BLOCKSIZE>
__global__ void compute_dh(double *d_dh, double *d_u, double *d_v, double H, double dx, double dy, int nx, int ny) {
    const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
    const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
    if (i < nx - 1 && j < ny - 1) { // avoid out-of-bounds access
        double du_dx = (d_u[(i + 1) * ny + j] - d_u[i * ny + j]) / dx;
        double dv_dy = (d_v[i * (ny + 1) + j + 1] - d_v[i * (ny + 1) + j]) / dy;
        d_dh[i * ny + j] = -H * (du_dx + dv_dy);
    }
}

// kernel to compute du on gpu w boundary checks
template <const uint BLOCKSIZE>
__global__ void compute_du(double *d_du, double *d_h, double g, double dx, int nx, int ny) {
    const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
    const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
    if (i < nx - 1 && j < ny) { // avoid out-of-bounds access
        double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - d_h[i * (ny + 1) + j]) / dx;
        d_du[i * ny + j] = -g * dh_dx;
    }
}

// kernel to compute dv on gpu w boundary checks
template <const uint BLOCKSIZE>
__global__ void compute_dv(double *d_dv, double *d_h, double g, double dy, int nx, int ny) {
    const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
    const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
    if (i < nx && j < ny - 1) { // avoid out-of-bounds access
        double dh_dy = (d_h[i * (ny + 1) + j + 1] - d_h[i * (ny + 1) + j]) / dy;
        d_dv[i * ny + j] = -g * dh_dy;
    }
}

// kernel to perform multistep update for h, u, v w additional boundary checks
template <const uint BLOCKSIZE>
__global__ void multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
                          double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
                          double a1, double a2, double a3, double dt, int nx, int ny) {
    const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
    const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
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

// Kernel to compute horizontal ghost cells
__global__ void compute_ghost_horizontal(double *d_h, int nx, int ny) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j < ny) {
        d_h[nx * (ny + 1) + j] = d_h[j];  // Copy the first row to the ghost row at the end
    }
}

// Kernel to compute vertical ghost cells
__global__ void compute_ghost_vertical(double *d_h, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nx) {
        d_h[i * (ny + 1) + ny] = d_h[i * (ny + 1)];  // Copy the first column to the ghost column at the end
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
    dim3 gridDim(ceil(nx / 32), ceil(ny / 32));
    dim3 blockDim(32 * 32);

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

    compute_ghost_horizontal<<<(ny + 31) / 32, 32>>>(d_h, nx, ny);
    compute_ghost_vertical<<<(nx + 31) / 32, 32>>>(d_h, nx, ny);
    cudaDeviceSynchronize();

    // launch kernels to compute derivs
    compute_dh<32><<<gridDim, blockDim>>>(d_dh, d_u, d_v, H, dx, dy, nx, ny);
    
    compute_du<32><<<gridDim, blockDim>>>(d_du, d_h, g, dx, nx, ny);

    compute_dv<32><<<gridDim, blockDim>>>(d_dv, d_h, g, dy, nx, ny);

    // launch kernel to update h, u, v with multistep method
    multistep<32><<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2, a1, a2, a3, dt, nx, ny);

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