// start original
// Timing gpu size 100
// Total time: 1.3071490000000001
// Timing gpu size 500
// Total time: 1.1688779999999999
// Timing gpu size 1000
// Total time: 2.319379
// Timing gpu size 2000
// Total time: 6.7092860000000005
// Timing gpu size 3000
// Total time: 13.767254
// Timing gpu size 4000
// Total time: 23.082241
// Timing gpu size 5000
// Total time: 37.448779
// Timing gpu size 6000
// Total time: 51.499084
// Timing gpu size 7000
// Total time: 69.772802
// Timing gpu size 8000
// Total time: 88.91878200000001
// Timing gpu size 9000
// Total time: 115.70036900000001
// Timing gpu size 10000
// Total time: 141.711399
// #include <cuda.h>
// #include <cuda_runtime.h>
// #include <math.h>
// #include <cstdio>     // fprintf
// #include <cstdlib>    // exit
// #include "../common/common.hpp"
// #include "../common/solver.hpp"
// #include <cublas_v2.h>

// // vars for grid size
// int nx, ny;
// double H, g, dx, dy, dt;

// // dev ptrs for fields and derivs
// double *d_h, *d_u, *d_v, *d_dh, *d_du, *d_dv, *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
// int t = 0;

// // kernel to compute dh on gpu w boundary checks
// template <const uint BLOCKSIZE>
// __global__ void compute_dh(double *d_dh, double *d_u, double *d_v, double H, double dx, double dy, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
//     if (i < nx - 1 && j < ny - 1) { // avoid out-of-bounds access
//         double du_dx = (d_u[(i + 1) * ny + j] - d_u[i * ny + j]) / dx;
//         double dv_dy = (d_v[i * (ny + 1) + j + 1] - d_v[i * (ny + 1) + j]) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }
// }

// // kernel to compute du on gpu w boundary checks
// template <const uint BLOCKSIZE>
// __global__ void compute_du(double *d_du, double *d_h, double g, double dx, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
//     if (i < nx - 1 && j < ny) { // avoid out-of-bounds access
//         double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - d_h[i * (ny + 1) + j]) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }
// }

// // kernel to compute dv on gpu w boundary checks
// template <const uint BLOCKSIZE>
// __global__ void compute_dv(double *d_dv, double *d_h, double g, double dy, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
//     if (i < nx && j < ny - 1) { // avoid out-of-bounds access
//         double dh_dy = (d_h[i * (ny + 1) + j + 1] - d_h[i * (ny + 1) + j]) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
// }

// // kernel to perform multistep update for h, u, v w additional boundary checks
// template <const uint BLOCKSIZE>
// __global__ void multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                           double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                           double a1, double a2, double a3, double dt, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
//     if (i < nx && j < ny) { // Ensure within bounds
//         d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
        
//         // Add boundary checks for d_u and d_v accesses
//         if (i + 1 < nx) { // Avoid out-of-bounds access for d_u
//             d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//         }
//         if (j + 1 < ny) { // Avoid out-of-bounds access for d_v
//             d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//         }
//     }
// }

// // kernel to compute horizontal ghost cells
// __global__ void compute_ghost_horizontal(double *d_h, int nx, int ny) {
//     int j = blockIdx.x * blockDim.x + threadIdx.x;
//     if (j < ny) {
//         d_h[nx * (ny + 1) + j] = d_h[j];  // Copy the first row to the ghost row at the end
//     }
// }

// // kernel to compute vertical ghost cells
// __global__ void compute_ghost_vertical(double *d_h, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     if (i < nx) {
//         d_h[i * (ny + 1) + ny] = d_h[i * (ny + 1)];  // Copy the first column to the ghost column at the end
//     }
// }

// // init fn to allocate mem on gpu and cp initial data 
// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     // set grid parameters
//     nx = nx_;
//     ny = ny_;
//     H = H_;
//     g = g_;
//     dx = length_ / nx;
//     dy = width_ / ny;
//     dt = dt_;

//     // alloc device mem
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

//     // cp init data to gpu
//     cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_u, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_v, v0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
// }

// // step fn to perform single time step on gpu 
// void step()
// {
//     // set up thread block and grid dimensions
//     dim3 gridDim(ceil(nx / 32), ceil(ny / 32));
//     dim3 blockDim(32 * 32);

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

//     compute_ghost_horizontal<<<(ny + 31) / 32, 32>>>(d_h, nx, ny);
//     compute_ghost_vertical<<<(nx + 31) / 32, 32>>>(d_h, nx, ny);
//     cudaDeviceSynchronize();

//     // launch kernels to compute derivs
//     compute_dh<32><<<gridDim, blockDim>>>(d_dh, d_u, d_v, H, dx, dy, nx, ny);
    
//     compute_du<32><<<gridDim, blockDim>>>(d_du, d_h, g, dx, nx, ny);

//     compute_dv<32><<<gridDim, blockDim>>>(d_dv, d_h, g, dy, nx, ny);

//     // launch kernel to update h, u, v with multistep method
//     multistep<32><<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2, a1, a2, a3, dt, nx, ny);

//     // swap deriv buffers
//     double *tmp;
//     tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
//     tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
//     tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;

//     t++;
// }

// // transfer fn to copy h field back to host 
// void transfer(double *h_host)
// {
//     cudaMemcpy(h_host, d_h, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyDeviceToHost);
// }

// // free gpu mem 
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
// // done original


// faster but fucked error
// Timing gpu size 100
// Total time: 1.2491999999999999
// Timing gpu size 500
// Total time: 1.159062
// Timing gpu size 1000
// Total time: 1.957747
// Timing gpu size 2000
// Total time: 5.112661999999999
// Timing gpu size 3000
// Total time: 10.386254
// Timing gpu size 4000
// Total time: 16.603002
// Timing gpu size 5000
// Total time: 26.914241
// Timing gpu size 6000
// Total time: 36.924604
// Timing gpu size 7000
// Total time: 50.473049
// Timing gpu size 8000
// Total time: 63.891193
// Timing gpu size 9000
// Total time: 84.134928
// Timing gpu size 10000
// Total time: 101.8472
// #include <cuda.h>
// #include <cuda_runtime.h>
// #include <math.h>
// #include <cstdio>
// #include <cstdlib>
// #include "../common/common.hpp"
// #include "../common/solver.hpp"
// #include <cublas_v2.h>

// // vars for grid size
// int nx, ny;
// double H, g, dx, dy, dt;

// // dev ptrs for fields and derivs
// double *d_h, *d_u, *d_v, *d_dh, *d_du, *d_dv, *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
// int t = 0;

// #define CUDA_CHECK(call) \
//     do { \
//         cudaError_t err = call; \
//         if (err != cudaSuccess) { \
//             fprintf(stderr, "CUDA error in %s at line %d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
//             exit(err); \
//         } \
//     } while (0)

// template <const uint BLOCKSIZE>
// __global__ void compute_and_update(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                                    double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                                    double H, double g, double dx, double dy, double a1, double a2, double a3, double dt, int nx, int ny) {
//     // Shared memory for local tiles
//     __shared__ double sh_h[BLOCKSIZE + 1][BLOCKSIZE + 1];
//     __shared__ double sh_u[BLOCKSIZE][BLOCKSIZE];
//     __shared__ double sh_v[BLOCKSIZE][BLOCKSIZE + 1];

//     // Calculate global and block indices
//     int i = blockIdx.x * BLOCKSIZE + threadIdx.y;
//     int j = blockIdx.y * BLOCKSIZE + threadIdx.x;

//     // Load data into shared memory with boundary padding
//     if (i < nx && j < ny) {
//         sh_h[threadIdx.y][threadIdx.x] = d_h[i * (ny + 1) + j];
//         if (threadIdx.x < BLOCKSIZE) sh_u[threadIdx.y][threadIdx.x] = d_u[i * ny + j];
//         if (threadIdx.y < BLOCKSIZE) sh_v[threadIdx.y][threadIdx.x] = d_v[i * (ny + 1) + j];

//         // Load boundary cells for shared memory padding
//         if (threadIdx.x == BLOCKSIZE - 1 && j + 1 < ny) sh_h[threadIdx.y][BLOCKSIZE] = d_h[i * (ny + 1) + j + 1];
//         if (threadIdx.y == BLOCKSIZE - 1 && i + 1 < nx) sh_h[BLOCKSIZE][threadIdx.x] = d_h[(i + 1) * (ny + 1) + j];

//         __syncthreads();

//         // Compute derivatives and update values in a single pass
//         if (i < nx - 1 && j < ny - 1) {
//             double du_dx = (sh_u[threadIdx.y][threadIdx.x + 1] - sh_u[threadIdx.y][threadIdx.x]) / dx;
//             double dv_dy = (sh_v[threadIdx.y + 1][threadIdx.x] - sh_v[threadIdx.y][threadIdx.x]) / dy;
//             d_dh[i * ny + j] = -H * (du_dx + dv_dy);

//             if (j < ny) {
//                 double dh_dx = (sh_h[threadIdx.y + 1][threadIdx.x] - sh_h[threadIdx.y][threadIdx.x]) / dx;
//                 d_du[i * ny + j] = -g * dh_dx;
//             }

//             if (i < nx) {
//                 double dh_dy = (sh_h[threadIdx.y][threadIdx.x + 1] - sh_h[threadIdx.y][threadIdx.x]) / dy;
//                 d_dv[i * ny + j] = -g * dh_dy;
//             }
//         }
//         __syncthreads();

//         // Update values using the multistep method
//         if (i < nx && j < ny) {
//             d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;

//             if (i + 1 < nx) {
//                 d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//             }
//             if (j + 1 < ny) {
//                 d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//             }
//         }
//     }
// }

// #define BLOCKSIZE 32

// void step()
// {
//     dim3 gridDim((nx + BLOCKSIZE - 1) / BLOCKSIZE, (ny + BLOCKSIZE - 1) / BLOCKSIZE);
//     dim3 blockDim(BLOCKSIZE, BLOCKSIZE);

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

//     compute_and_update<32><<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2,
//                                                          H, g, dx, dy, a1, a2, a3, dt, nx, ny);
//     cudaDeviceSynchronize();

//     // Swap derivative buffers
//     double *tmp;
//     tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
//     tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
//     tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;

//     t++;
// }

// // transfer fn to copy h field back to host 
// void transfer(double *h_host)
// {
//     cudaMemcpy(h_host, d_h, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyDeviceToHost);
// }

// // init fn to allocate mem on gpu and cp initial data 
// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     // set grid parameters
//     nx = nx_;
//     ny = ny_;
//     H = H_;
//     g = g_;
//     dx = length_ / nx;
//     dy = width_ / ny;
//     dt = dt_;

//     // alloc device mem
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

//     // cp init data to gpu
//     cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_u, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_v, v0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
// }

// // free gpu mem 
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





// before combining kernels
// #include <cuda.h>
// #include <cuda_runtime.h>
// #include <math.h>
// #include <cstdio>     // fprintf
// #include <cstdlib>    // exit
// #include "../common/common.hpp"
// #include "../common/solver.hpp"
// #include <cublas_v2.h>

// // vars for grid size
// int nx, ny;
// double H, g, dx, dy, dt;

// // dev ptrs for fields and derivs
// double *d_h, *d_u, *d_v, *d_dh, *d_du, *d_dv, *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
// int t = 0;

// // kernel to compute dh on gpu w boundary checks
// template <const uint BLOCKSIZE>
// __global__ void compute_dh(double *d_dh, double *d_u, double *d_v, double H, double dx, double dy, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
//     if (i < nx - 1 && j < ny - 1) { // avoid out-of-bounds access
//         double du_dx = (d_u[(i + 1) * ny + j] - d_u[i * ny + j]) / dx;
//         double dv_dy = (d_v[i * (ny + 1) + j + 1] - d_v[i * (ny + 1) + j]) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }
// }

// // kernel to compute du on gpu w boundary checks
// template <const uint BLOCKSIZE>
// __global__ void compute_du(double *d_du, double *d_h, double g, double dx, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
//     if (i < nx - 1 && j < ny) { // avoid out-of-bounds access
//         double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - d_h[i * (ny + 1) + j]) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }
// }

// // kernel to compute dv on gpu w boundary checks
// template <const uint BLOCKSIZE>
// __global__ void compute_dv(double *d_dv, double *d_h, double g, double dy, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
//     if (i < nx && j < ny - 1) { // avoid out-of-bounds access
//         double dh_dy = (d_h[i * (ny + 1) + j + 1] - d_h[i * (ny + 1) + j]) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
// }

// // kernel to perform multistep update for h, u, v w additional boundary checks
// template <const uint BLOCKSIZE>
// __global__ void multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                           double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                           double a1, double a2, double a3, double dt, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
//     if (i < nx && j < ny) { // Ensure within bounds
//         d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
        
//         // Add boundary checks for d_u and d_v accesses
//         if (i + 1 < nx) { // Avoid out-of-bounds access for d_u
//             d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//         }
//         if (j + 1 < ny) { // Avoid out-of-bounds access for d_v
//             d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//         }
//     }
// }

// // kernel to compute horizontal ghost cells
// __global__ void compute_ghost_horizontal(double *d_h, int nx, int ny) {
//     int j = blockIdx.x * blockDim.x + threadIdx.x;
//     if (j < ny) {
//         d_h[nx * (ny + 1) + j] = d_h[j];  // Copy the first row to the ghost row at the end
//     }
// }

// // kernel to compute vertical ghost cells
// __global__ void compute_ghost_vertical(double *d_h, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     if (i < nx) {
//         d_h[i * (ny + 1) + ny] = d_h[i * (ny + 1)];  // Copy the first column to the ghost column at the end
//     }
// }

// // kernel to compute boundaries for horizontal velocity u
// __global__ void compute_boundaries_horizontal(double *d_u, int nx, int ny) {
//     int j = blockIdx.x * blockDim.x + threadIdx.x;
//     if (j < ny) {
//         d_u[0 * ny + j] = d_u[(nx - 1) * ny + j];  // Wrap boundary
//     }
// }

// // kernel to compute boundaries for vertical velocity v
// __global__ void compute_boundaries_vertical(double *d_v, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     if (i < nx) {
//         d_v[i * (ny + 1) + 0] = d_v[i * (ny + 1) + (ny - 1)];  // Wrap boundary
//     }
// }

// // init fn to allocate mem on gpu and cp initial data 
// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     // set grid parameters
//     nx = nx_;
//     ny = ny_;
//     H = H_;
//     g = g_;
//     dx = length_ / nx;
//     dy = width_ / ny;
//     dt = dt_;

//     // alloc device mem
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

//     // cp init data to gpu
//     cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_u, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_v, v0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
// }

// // step fn to perform single time step on gpu 
// void step()
// {
//     // set up thread block and grid dimensions
//     dim3 gridDim(ceil(nx / 32), ceil(ny / 32));
//     dim3 blockDim(32 * 32);

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

//     compute_ghost_horizontal<<<(ny + 31) / 32, 32>>>(d_h, nx, ny);
//     compute_ghost_vertical<<<(nx + 31) / 32, 32>>>(d_h, nx, ny);

//     compute_boundaries_horizontal<<<(ny + 31) / 32, 32>>>(d_u, nx, ny);
//     compute_boundaries_vertical<<<(nx + 31) / 32, 32>>>(d_v, nx, ny);
    
//     cudaDeviceSynchronize();

//     // launch kernels to compute derivs
//     compute_dh<32><<<gridDim, blockDim>>>(d_dh, d_u, d_v, H, dx, dy, nx, ny);
    
//     compute_du<32><<<gridDim, blockDim>>>(d_du, d_h, g, dx, nx, ny);

//     compute_dv<32><<<gridDim, blockDim>>>(d_dv, d_h, g, dy, nx, ny);

//     // launch kernel to update h, u, v with multistep method
//     multistep<32><<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2, a1, a2, a3, dt, nx, ny);

//     // swap deriv buffers
//     double *tmp;
//     tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
//     tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
//     tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;

//     t++;
// }

// // transfer fn to copy h field back to host 
// void transfer(double *h_host)
// {
//     cudaMemcpy(h_host, d_h, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyDeviceToHost);
// }

// // free gpu mem 
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







// before shared mem
// #include <cuda.h>
// #include <cuda_runtime.h>
// #include <math.h>
// #include <cstdio>
// #include <cstdlib>
// #include "../common/common.hpp"
// #include "../common/solver.hpp"
// #include <cublas_v2.h>

// // vars for grid size
// int nx, ny;
// double H, g, dx, dy, dt;

// // dev ptrs for fields and derivs
// double *d_h, *d_u, *d_v, *d_dh, *d_du, *d_dv, *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
// int t = 0;

// // Device function to compute ghost cells for h
// __device__ void compute_ghost_cells(double *d_h, int nx, int ny, int j, int i) {
//     if (j < ny && i == nx) {
//         d_h[i * (ny + 1) + j] = d_h[j];  // Horizontal ghost cell at the boundary
//     }
//     if (i < nx && j == ny) {
//         d_h[i * (ny + 1) + j] = d_h[i * (ny + 1)];  // Vertical ghost cell at the boundary
//     }
// }

// // Device function to apply boundary conditions for u and v
// __device__ void apply_boundaries(double *d_u, double *d_v, int nx, int ny, int j, int i) {
//     if (j < ny && i == 0) {
//         d_u[i * ny + j] = d_u[(nx - 1) * ny + j];  // Horizontal boundary for u
//     }
//     if (i < nx && j == 0) {
//         d_v[i * (ny + 1) + j] = d_v[i * (ny + 1) + (ny - 1)];  // Vertical boundary for v
//     }
// }

// // Device function to compute dh
// __device__ void compute_dh(double *d_dh, double *d_u, double *d_v, double H, double dx, double dy, int nx, int ny, int j, int i) {
//     if (i < nx - 1 && j < ny - 1) {
//         double du_dx = (d_u[(i + 1) * ny + j] - d_u[i * ny + j]) / dx;
//         double dv_dy = (d_v[i * (ny + 1) + j + 1] - d_v[i * (ny + 1) + j]) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }
// }

// // Device function to compute du
// __device__ void compute_du(double *d_du, double *d_h, double g, double dx, int nx, int ny, int j, int i) {
//     if (i < nx - 1 && j < ny) {
//         double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - d_h[i * (ny + 1) + j]) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }
// }

// // Device function to compute dv
// __device__ void compute_dv(double *d_dv, double *d_h, double g, double dy, int nx, int ny, int j, int i) {
//     if (i < nx && j < ny - 1) {
//         double dh_dy = (d_h[i * (ny + 1) + j + 1] - d_h[i * (ny + 1) + j]) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
// }

// // Combined kernel to compute ghost cells, boundaries, and derivatives
// template <const uint BLOCKSIZE>
// __global__ void compute_step(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                              double H, double g, double dx, double dy, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);

//     // Call individual functions to perform tasks
//     compute_ghost_cells(d_h, nx, ny, j, i);
//     apply_boundaries(d_u, d_v, nx, ny, j, i);
//     compute_dh(d_dh, d_u, d_v, H, dx, dy, nx, ny, j, i);
//     compute_du(d_du, d_h, g, dx, nx, ny, j, i);
//     compute_dv(d_dv, d_h, g, dy, nx, ny, j, i);
// }

// // Kernel to perform multistep update for h, u, v
// template <const uint BLOCKSIZE>
// __global__ void multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                           double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                           double a1, double a2, double a3, double dt, int nx, int ny) {
//     const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
//     const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
//     if (i < nx && j < ny) {
//         d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
//         if (i + 1 < nx) {
//             d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//         }
//         if (j + 1 < ny) {
//             d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//         }
//     }
// }

// // Initialize GPU memory and copy initial data
// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_) {
//     nx = nx_;
//     ny = ny_;
//     H = H_;
//     g = g_;
//     dx = length_ / nx;
//     dy = width_ / ny;
//     dt = dt_;

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

//     cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_u, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_v, v0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
// }

// // Step function to perform a single time step on the GPU
// void step() {
//     dim3 gridDim(ceil(nx / 32), ceil(ny / 32));
//     dim3 blockDim(32 * 32);

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

//     compute_step<32><<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, H, g, dx, dy, nx, ny);
//     cudaDeviceSynchronize();

//     multistep<32><<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2, a1, a2, a3, dt, nx, ny);
//     cudaDeviceSynchronize();

//     double *tmp;
//     tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
//     tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
//     tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;

//     t++;
// }

// // Transfer function to copy h field back to host
// void transfer(double *h_host) {
//     cudaMemcpy(h_host, d_h, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyDeviceToHost);
// }

// // Free GPU memory
// void free_memory() {
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



// Timing gpu size 100
// Total time: 1.2632009999999998
// Timing gpu size 500
// Total time: 1.3242990000000001
// Timing gpu size 1000
// Total time: 2.260982
// Timing gpu size 2000
// Total time: 5.831777
// Timing gpu size 3000
// Total time: 11.737845
// Timing gpu size 4000
// Total time: 20.241487000000003
// Timing gpu size 5000
// Total time: 31.231681
// Timing gpu size 6000
// Total time: 42.637890999999996
// Timing gpu size 7000
// Total time: 57.523322
// Timing gpu size 8000
// Total time: 76.918652
// Timing gpu size 9000
// Total time: 95.27791400000001
// Timing gpu size 10000
// Total time: 116.987207
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include "../common/common.hpp"
#include "../common/solver.hpp"
#include <cublas_v2.h>

// vars for grid size
int nx, ny;
double H, g, dx, dy, dt;

// dev ptrs for fields and derivs
double *d_h, *d_u, *d_v, *d_dh, *d_du, *d_dv, *d_dh1, *d_du1, *d_dv1, *d_dh2, *d_du2, *d_dv2;
int t = 0;

// Define TILE_SIZE based on Perlmutter's shared memory limit
#define TILE_SIZE 32

// Combined kernel with shared memory for efficient data access
template <const uint BLOCKSIZE>
__global__ void compute_step(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
                             double H, double g, double dx, double dy, int nx, int ny) {

    // Define shared memory for each field, with padding for boundary cells
    __shared__ double sh_h[TILE_SIZE + 1][TILE_SIZE + 1];
    __shared__ double sh_u[TILE_SIZE + 1][TILE_SIZE];
    __shared__ double sh_v[TILE_SIZE][TILE_SIZE + 1];

    // Calculate 2D indices within the block using a 1D thread layout
    const int local_i = threadIdx.x / BLOCKSIZE;   // Row index within tile
    const int local_j = threadIdx.x % BLOCKSIZE;   // Column index within tile

    // Calculate global indices
    const int i = blockIdx.x * TILE_SIZE + local_i;
    const int j = blockIdx.y * TILE_SIZE + local_j;

    // Load data into shared memory, with boundary checks
    if (i <= nx && j <= ny) {
        sh_h[local_i][local_j] = d_h[i * (ny + 1) + j];
    }
    if (i < nx && j < ny) {
        sh_u[local_i][local_j] = d_u[i * ny + j];
        sh_v[local_i][local_j] = d_v[i * (ny + 1) + j];
    }

    // Load extra data for ghost cells and boundaries
    if (local_i == TILE_SIZE - 1 && i + 1 <= nx) {
        sh_h[local_i + 1][local_j] = d_h[(i + 1) * (ny + 1) + j];
    }
    if (local_j == TILE_SIZE - 1 && j + 1 <= ny) {
        sh_h[local_i][local_j + 1] = d_h[i * (ny + 1) + (j + 1)];
    }

    __syncthreads();

    // Compute ghost cells for h
    if (local_j < TILE_SIZE && local_i == TILE_SIZE - 1 && i == nx) {
        sh_h[local_i + 1][local_j] = sh_h[0][local_j];
    }
    if (local_i < TILE_SIZE && local_j == TILE_SIZE - 1 && j == ny) {
        sh_h[local_i][local_j + 1] = sh_h[local_i][0];
    }

    // Apply boundary conditions for u and v
    if (j < ny && i == 0) {
        sh_u[local_i][local_j] = sh_u[TILE_SIZE - 1][local_j];  // Horizontal boundary for u
    }
    if (i < nx && j == 0) {
        sh_v[local_i][local_j] = sh_v[local_i][TILE_SIZE - 1];  // Vertical boundary for v
    }

    __syncthreads();

    // Compute dh using shared memory
    if (i < nx - 1 && j < ny - 1) {
        double du_dx = (sh_u[local_i + 1][local_j] - sh_u[local_i][local_j]) / dx;
        double dv_dy = (sh_v[local_i][local_j + 1] - sh_v[local_i][local_j]) / dy;
        d_dh[i * ny + j] = -H * (du_dx + dv_dy);
    }

    // Compute du using shared memory
    if (i < nx - 1 && j < ny) {
        double dh_dx = (sh_h[local_i + 1][local_j] - sh_h[local_i][local_j]) / dx;
        d_du[i * ny + j] = -g * dh_dx;
    }

    // Compute dv using shared memory
    if (i < nx && j < ny - 1) {
        double dh_dy = (sh_h[local_i][local_j + 1] - sh_h[local_i][local_j]) / dy;
        d_dv[i * ny + j] = -g * dh_dy;
    }
}

// Kernel to perform multistep update for h, u, v
template <const uint BLOCKSIZE>
__global__ void multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
                          double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
                          double a1, double a2, double a3, double dt, int nx, int ny) {
    const int i = blockIdx.x * BLOCKSIZE + (threadIdx.x / BLOCKSIZE);
    const int j = blockIdx.y * BLOCKSIZE + (threadIdx.x % BLOCKSIZE);
    
    if (i < nx && j < ny) {
        d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
        if (i + 1 < nx) {
            d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
        }
        if (j + 1 < ny) {
            d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
        }
    }
}

// Initialize GPU memory and copy initial data
void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_) {
    nx = nx_;
    ny = ny_;
    H = H_;
    g = g_;
    dx = length_ / nx;
    dy = width_ / ny;
    dt = dt_;

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

    cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_u, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, v0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
}

// Step function to perform a single time step on the GPU
void step() {
    dim3 gridDim(ceil(nx / TILE_SIZE), ceil(ny / TILE_SIZE));
    dim3 blockDim(TILE_SIZE, TILE_SIZE);

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

    compute_step<32><<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, H, g, dx, dy, nx, ny);
    cudaDeviceSynchronize();

    multistep<32><<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2, a1, a2, a3, dt, nx, ny);
    cudaDeviceSynchronize();

    double *tmp;
    tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
    tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
    tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;

    t++;
}

// Transfer function to copy h field back to host
void transfer(double *h_host) {
    cudaMemcpy(h_host, d_h, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyDeviceToHost);
}

// Free GPU memory
void free_memory() {
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
