//Max error: 2.220446049250313e-14
// Timing gpu size 100
// Total time: 1.183839
// Timing gpu size 500
// Total time: 1.0986870000000002
// Timing gpu size 1000
// Total time: 2.277707
// Timing gpu size 2000
// Total time: 5.895163999999999
// Timing gpu size 3000
// Total time: 12.073941999999999
// Timing gpu size 4000
// Total time: 20.583716000000003
// Timing gpu size 5000
// Total time: 31.6028
// Timing gpu size 6000
// Total time: 44.914028
// Timing gpu size 7000
// Total time: 60.912118
// Timing gpu size 8000
// Total time: 78.97792899999999
// Timing gpu size 9000
// Total time: 99.930127
// Timing gpu size 10000
// Total time: 122.842415
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

// // Combined kernel to compute ghost cells, boundaries, and derivatives
// __global__ void compute_step(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                              double H, double g, double dx, double dy, int nx, int ny) {
//     int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     int i = idx / ny; // Convert linear idx to 2D grid coordinates
//     int j = idx % ny;

//     if (i >= nx || j >= ny) return; // Boundary check

//     // Calculate ghost cells for h
//     if (j < ny && i == nx) d_h[i * (ny + 1) + j] = d_h[j];
//     if (i < nx && j == ny) d_h[i * (ny + 1) + j] = d_h[i * (ny + 1)];

//     // Apply boundary conditions for u and v
//     if (j < ny && i == 0) d_u[i * ny + j] = d_u[(nx - 1) * ny + j];
//     if (i < nx && j == 0) d_v[i * (ny + 1) + j] = d_v[i * (ny + 1) + (ny - 1)];

//     // Compute dh
//     if (i < nx - 1 && j < ny - 1) {
//         double du_dx = (d_u[(i + 1) * ny + j] - d_u[i * ny + j]) / dx;
//         double dv_dy = (d_v[i * (ny + 1) + j + 1] - d_v[i * (ny + 1) + j]) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }

//     // Compute du
//     if (i < nx - 1 && j < ny) {
//         double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - d_h[i * (ny + 1) + j]) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }

//     // Compute dv
//     if (i < nx && j < ny - 1) {
//         double dh_dy = (d_h[i * (ny + 1) + j + 1] - d_h[i * (ny + 1) + j]) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
// }

// // Kernel to perform multistep update for h, u, v
// __global__ void multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                           double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                           double a1, double a2, double a3, double dt, int nx, int ny) {
//     int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     int i = idx / ny;
//     int j = idx % ny;

//     if (i >= nx || j >= ny) return;

//     d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
//     if (i + 1 < nx) {
//         d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//     }
//     if (j + 1 < ny) {
//         d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//     }
// }

// // Step function to perform a single time step on the GPU
// void step() {
//     dim3 gridDim((nx * ny + 1023) / 1024); // Launch enough blocks to cover all elements
//     dim3 blockDim(1024); // 1024 threads per block

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

//     compute_step<<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, H, g, dx, dy, nx, ny);
//     cudaDeviceSynchronize();

//     multistep<<<gridDim, blockDim>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2, a1, a2, a3, dt, nx, ny);
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


// wip
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

// #define BLOCK_DIM_X 16  // Block width
// #define BLOCK_DIM_Y 16  // Block height
// #define BLOCKSIZE (BLOCK_DIM_X * BLOCK_DIM_Y)

// #define BLOCK_DIM_X 16  // Block width
// #define BLOCK_DIM_Y 16  // Block height
// #define BLOCKSIZE (BLOCK_DIM_X * BLOCK_DIM_Y)

// __global__ void compute_and_multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                                       double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                                       double H, double g, double dx, double dy, double a1, double a2, double a3, double dt,
//                                       int nx, int ny) {
//     extern __shared__ double shared_mem[];
//     double *h_shared = shared_mem;
//     double *u_shared = h_shared + (BLOCK_DIM_X + 2) * (BLOCK_DIM_Y + 2);
//     double *v_shared = u_shared + (BLOCK_DIM_X + 2) * (BLOCK_DIM_Y + 2);

//     int tx = threadIdx.x;
//     int ty = threadIdx.y;
//     int i = blockIdx.x * BLOCK_DIM_X + tx;  // Global x index
//     int j = blockIdx.y * BLOCK_DIM_Y + ty;  // Global y index

//     // Shared memory indexing with padding (1-cell ghost padding on each side)
//     int shared_idx = (ty + 1) * (BLOCK_DIM_X + 2) + (tx + 1);

//     // Load the main cell and boundary cells into shared memory with padding
//     if (i < nx && j < ny) {
//         h_shared[shared_idx] = d_h[i * (ny + 1) + j];
//         u_shared[shared_idx] = d_u[i * ny + j];
//         v_shared[shared_idx] = d_v[i * (ny + 1) + j];
//     }

//     // Load ghost cells along each side
//     if (tx == 0 && i > 0) {
//         h_shared[shared_idx - 1] = d_h[(i - 1) * (ny + 1) + j];
//         u_shared[shared_idx - 1] = d_u[(i - 1) * ny + j];
//         v_shared[shared_idx - 1] = d_v[(i - 1) * (ny + 1) + j];
//     }
//     if (tx == BLOCK_DIM_X - 1 && i < nx - 1) {
//         h_shared[shared_idx + 1] = d_h[(i + 1) * (ny + 1) + j];
//         u_shared[shared_idx + 1] = d_u[(i + 1) * ny + j];
//         v_shared[shared_idx + 1] = d_v[(i + 1) * (ny + 1) + j];
//     }
//     if (ty == 0 && j > 0) {
//         h_shared[shared_idx - (BLOCK_DIM_X + 2)] = d_h[i * (ny + 1) + j - 1];
//         u_shared[shared_idx - (BLOCK_DIM_X + 2)] = d_u[i * ny + j - 1];
//         v_shared[shared_idx - (BLOCK_DIM_X + 2)] = d_v[i * (ny + 1) + j - 1];
//     }
//     if (ty == BLOCK_DIM_Y - 1 && j < ny - 1) {
//         h_shared[shared_idx + (BLOCK_DIM_X + 2)] = d_h[i * (ny + 1) + j + 1];
//         u_shared[shared_idx + (BLOCK_DIM_X + 2)] = d_u[i * ny + j + 1];
//         v_shared[shared_idx + (BLOCK_DIM_X + 2)] = d_v[i * (ny + 1) + j + 1];
//     }
//     __syncthreads();

//     // Global boundary check
//     if (i >= nx || j >= ny) return;

//     // Compute dh: finite differences along x and y using shared memory
//     if (tx < BLOCK_DIM_X && ty < BLOCK_DIM_Y) {
//         double du_dx = (u_shared[shared_idx + (BLOCK_DIM_X + 2)] - u_shared[shared_idx]) / dx;
//         double dv_dy = (v_shared[shared_idx + 1] - v_shared[shared_idx]) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }

//     // Compute du: finite differences along x using shared memory
//     if (tx < BLOCK_DIM_X && ty < BLOCK_DIM_Y) {
//         double dh_dx = (h_shared[shared_idx + (BLOCK_DIM_X + 2)] - h_shared[shared_idx]) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }

//     // Compute dv: finite differences along y using shared memory
//     if (tx < BLOCK_DIM_X && ty < BLOCK_DIM_Y) {
//         double dh_dy = (h_shared[shared_idx + 1] - h_shared[shared_idx]) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
//     __syncthreads();

//     // Multistep update for h, u, v
//     d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
//     if (i + 1 < nx) {
//         d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//     }
//     if (j + 1 < ny) {
//         d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//     }
// }

// void step() {
//     dim3 blockDim(BLOCK_DIM_X, BLOCK_DIM_Y);
//     dim3 gridDim((nx + BLOCK_DIM_X - 1) / BLOCK_DIM_X, (ny + BLOCK_DIM_Y - 1) / BLOCK_DIM_Y);

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

//     // int shared_mem_size = 3 * (BLOCK_DIM_X + 1) * (BLOCK_DIM_Y + 1) * sizeof(double);
//     int shared_mem_size = 3 * (BLOCK_DIM_X + 2) * (BLOCK_DIM_Y + 2) * sizeof(double);
//     compute_and_multistep<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv,
//                                                                   d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2,
//                                                                   H, g, dx, dy, a1, a2, a3, dt, nx, ny);
//     cudaDeviceSynchronize();

//     // Rotate pointers for multistep updates
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






//keep BEST!
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

// #define BLOCKSIZE 512  // Adjusted for better occupancy on Perlmutter

// WORKS, for 1D blocking
// __global__ void compute_and_multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                                       double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                                       double H, double g, double dx, double dy, double a1, double a2, double a3, double dt,
//                                       int nx, int ny) {
//     int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     int i = idx / ny; // Convert linear idx to 2D grid coordinates
//     int j = idx % ny;

//     if (i >= nx || j >= ny) return; // Boundary check

//     // Directly load values from global memory
//     double h = d_h[i * (ny + 1) + j];
//     double u = d_u[i * ny + j];
//     double v = d_v[i * (ny + 1) + j];

//     // Calculate ghost cells for h
//     if (j < ny && i == nx) d_h[i * (ny + 1) + j] = d_h[(nx - 1) * (ny + 1) + j];  // Last row
//     if (i < nx && j == ny) d_h[i * (ny + 1) + j] = d_h[i * (ny + 1) + (ny - 1)];  // Last column

//     // Apply boundary conditions for u and v
//     if (j < ny && i == 0) d_u[i * ny + j] = d_u[(nx - 1) * ny + j];  // First row
//     if (i < nx && j == 0) d_v[i * (ny + 1) + j] = d_v[i * (ny + 1) + (ny - 1)]; // First column

//     // Compute dh: finite differences along x and y using direct global memory access
//     if (i < nx - 1 && j < ny - 1) {
//         double du_dx = (d_u[(i + 1) * ny + j] - u) / dx;
//         double dv_dy = (d_v[i * (ny + 1) + j + 1] - v) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }

//     // Compute du: finite differences along x using direct global memory access
//     if (i < nx - 1 && j < ny) {
//         double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - h) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }

//     // Compute dv: finite differences along y using direct global memory access
//     if (i < nx && j < ny - 1) {
//         double dh_dy = (d_h[i * (ny + 1) + j + 1] - h) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }

//     __syncthreads();

//     // Multistep update for h, u, v
//     d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
//     if (i + 1 < nx) {
//         d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//     }
//     if (j + 1 < ny) {
//         d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//     }
// }

// // 1D
// void step() {
//     dim3 gridDim((nx * ny + BLOCKSIZE - 1) / BLOCKSIZE); // Launch enough blocks to cover all elements
//     dim3 blockDim(BLOCKSIZE);

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

//     int shared_mem_size = 3 * BLOCKSIZE * sizeof(double);
//     compute_and_multistep<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv,
//                                                                   d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2,
//                                                                   H, g, dx, dy, a1, a2, a3, dt, nx, ny);
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




//works also
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

// #define BLOCK_X 16
// #define BLOCK_Y 16

// __global__ void compute_derivatives(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                                     double H, double g, double dx, double dy, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;

//     if (i >= nx || j >= ny) return;

//     double h = d_h[i * (ny + 1) + j];
//     double u = d_u[i * ny + j];
//     double v = d_v[i * (ny + 1) + j];

//     // Calculate ghost cells for h
//     if (j < ny && i == nx) d_h[i * (ny + 1) + j] = d_h[(nx - 1) * (ny + 1) + j];
//     if (i < nx && j == ny) d_h[i * (ny + 1) + j] = d_h[i * (ny + 1) + (ny - 1)];

//     // Apply boundary conditions for u and v
//     if (j < ny && i == 0) d_u[i * ny + j] = d_u[(nx - 1) * ny + j];
//     if (i < nx && j == 0) d_v[i * (ny + 1) + j] = d_v[i * (ny + 1) + (ny - 1)];

//     // Compute derivatives as in the original kernel
//     if (i < nx - 1 && j < ny - 1) {
//         double du_dx = (d_u[(i + 1) * ny + j] - u) / dx;
//         double dv_dy = (d_v[i * (ny + 1) + j + 1] - v) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }

//     if (i < nx - 1 && j < ny) {
//         double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - h) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }

//     if (i < nx && j < ny - 1) {
//         double dh_dy = (d_h[i * (ny + 1) + j + 1] - h) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
// }

// __global__ void update_fields(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                               double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                               double a1, double a2, double a3, double dt, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;

//     if (i >= nx || j >= ny) return;

//     // Update h, u, and v as in the original kernel
//     d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
//     if (i + 1 < nx) {
//         d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//     }
//     if (j + 1 < ny) {
//         d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//     }
// }

// void step() {
//     dim3 blockDim(BLOCK_X, BLOCK_Y);
//     dim3 gridDim((nx + BLOCK_X - 1) / BLOCK_X, (ny + BLOCK_Y - 1) / BLOCK_Y);

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

//     int shared_mem_size = 3 * BLOCK_X * BLOCK_Y * sizeof(double);

//     // Launch kernel to compute derivatives
//     compute_derivatives<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, H, g, dx, dy, nx, ny);
//     cudaDeviceSynchronize();

//     // Launch kernel to update fields
//     update_fields<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv,
//                                                           d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2,
//                                                           a1, a2, a3, dt, nx, ny);
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



// keep, but slower
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

// #define BLOCK_X 16
// #define BLOCK_Y 16

// __global__ void compute_derivatives(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                                     double H, double g, double dx, double dy, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
//     int local_i = threadIdx.x + 1;  // Local index within shared memory, with padding for ghost cells
//     int local_j = threadIdx.y + 1;

//     extern __shared__ double shared_mem[];
//     double *sh_h = shared_mem;
//     double *sh_u = shared_mem + (BLOCK_X + 2) * (BLOCK_Y + 2);
//     double *sh_v = sh_u + (BLOCK_X + 2) * (BLOCK_Y + 2);

//     // Load h, u, v into shared memory with padding for ghost cells
//     if (i < nx && j < ny) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j] = d_h[i * (ny + 1) + j];
//         sh_u[local_i * (BLOCK_Y + 2) + local_j] = d_u[i * ny + j];
//         sh_v[local_i * (BLOCK_Y + 2) + local_j] = d_v[i * (ny + 1) + j];
//     }

//     // Load ghost cells for h, u, and v
//     if (threadIdx.x == 0 && i > 0) {
//         sh_h[(local_i - 1) * (BLOCK_Y + 2) + local_j] = d_h[(i - 1) * (ny + 1) + j];
//         sh_u[(local_i - 1) * (BLOCK_Y + 2) + local_j] = d_u[(i - 1) * ny + j];
//         sh_v[(local_i - 1) * (BLOCK_Y + 2) + local_j] = d_v[(i - 1) * (ny + 1) + j];
//     }
//     if (threadIdx.x == BLOCK_X - 1 && i < nx - 1) {
//         sh_h[(local_i + 1) * (BLOCK_Y + 2) + local_j] = d_h[(i + 1) * (ny + 1) + j];
//         sh_u[(local_i + 1) * (BLOCK_Y + 2) + local_j] = d_u[(i + 1) * ny + j];
//         sh_v[(local_i + 1) * (BLOCK_Y + 2) + local_j] = d_v[(i + 1) * (ny + 1) + j];
//     }
//     if (threadIdx.y == 0 && j > 0) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j - 1] = d_h[i * (ny + 1) + j - 1];
//         sh_u[local_i * (BLOCK_Y + 2) + local_j - 1] = d_u[i * ny + j - 1];
//         sh_v[local_i * (BLOCK_Y + 2) + local_j - 1] = d_v[i * (ny + 1) + j - 1];
//     }
//     if (threadIdx.y == BLOCK_Y - 1 && j < ny - 1) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j + 1] = d_h[i * (ny + 1) + j + 1];
//         sh_u[local_i * (BLOCK_Y + 2) + local_j + 1] = d_u[i * ny + j + 1];
//         sh_v[local_i * (BLOCK_Y + 2) + local_j + 1] = d_v[i * (ny + 1) + j + 1];
//     }

//     __syncthreads();

//     // Apply boundary conditions for h, u, v at the domain edges
//     if (i == nx) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j] = d_h[(nx - 1) * (ny + 1) + j];
//     }
//     if (j == ny) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j] = d_h[i * (ny + 1) + (ny - 1)];
//     }
//     if (i == 0) {
//         sh_u[local_i * BLOCK_Y + local_j] = d_u[(nx - 1) * ny + j];
//     }
//     if (j == 0) {
//         sh_v[local_i * (BLOCK_Y + 2) + local_j] = d_v[i * (ny + 1) + (ny - 1)];
//     }

//     __syncthreads();

//     // Compute derivatives using shared memory
//     if (i < nx - 1 && j < ny - 1) {
//         double du_dx = (sh_u[(local_i + 1) * (BLOCK_Y + 2) + local_j] - sh_u[local_i * (BLOCK_Y + 2) + local_j]) / dx;
//         double dv_dy = (sh_v[local_i * (BLOCK_Y + 2) + local_j + 1] - sh_v[local_i * (BLOCK_Y + 2) + local_j]) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }

//     if (i < nx - 1 && j < ny) {
//         double dh_dx = (sh_h[(local_i + 1) * (BLOCK_Y + 2) + local_j] - sh_h[local_i * (BLOCK_Y + 2) + local_j]) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }

//     if (i < nx && j < ny - 1) {
//         double dh_dy = (sh_h[local_i * (BLOCK_Y + 2) + local_j + 1] - sh_h[local_i * (BLOCK_Y + 2) + local_j]) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
// }


// __global__ void update_fields(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                               double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                               double a1, double a2, double a3, double dt, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;

//     if (i >= nx || j >= ny) return;

//     // Update h, u, and v as in the original kernel
//     d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
//     if (i + 1 < nx) {
//         d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//     }
//     if (j + 1 < ny) {
//         d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//     }
// }

// void step() {
//     dim3 blockDim(BLOCK_X, BLOCK_Y);
//     dim3 gridDim((nx + BLOCK_X - 1) / BLOCK_X, (ny + BLOCK_Y - 1) / BLOCK_Y);

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

//     int shared_mem_size = (BLOCK_X + 2) * (BLOCK_Y + 2) * 3 * sizeof(double); // Shared memory size for h, u, v

//     // Launch kernel to compute derivatives
//     compute_derivatives<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv, H, g, dx, dy, nx, ny);
//     cudaDeviceSynchronize();

//     // Launch kernel to update fields
//     update_fields<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv,
//                                                           d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2,
//                                                           a1, a2, a3, dt, nx, ny);
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





//the best NOW
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

#define BLOCK_X 16
#define BLOCK_Y 16

__global__ void compute_and_update(double *d_h, double *d_u, double *d_v,
                                   double *d_dh, double *d_du, double *d_dv,
                                   double *d_dh1, double *d_du1, double *d_dv1,
                                   double *d_dh2, double *d_du2, double *d_dv2,
                                   double H, double g, double dx, double dy,
                                   double a1, double a2, double a3, double dt,
                                   int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int local_i = threadIdx.x + 1;  // Local index within shared memory with padding for ghost cells
    int local_j = threadIdx.y + 1;

    extern __shared__ double shared_mem[];
    double *sh_h = shared_mem;
    double *sh_u = sh_h + (BLOCK_X + 2) * (BLOCK_Y + 2);
    double *sh_v = sh_u + (BLOCK_X + 2) * (BLOCK_Y + 2);

    // Load h, u, v into shared memory with padding for ghost cells
    if (i < nx && j < ny) {
        sh_h[local_i * (BLOCK_Y + 2) + local_j] = d_h[i * (ny + 1) + j];
        sh_u[local_i * (BLOCK_Y + 2) + local_j] = d_u[i * ny + j];
        sh_v[local_i * (BLOCK_Y + 2) + local_j] = d_v[i * (ny + 1) + j];
    }

    // Load halo regions for boundary values
    if (threadIdx.x == 0 && i > 0) {
        sh_h[(local_i - 1) * (BLOCK_Y + 2) + local_j] = d_h[(i - 1) * (ny + 1) + j];
        sh_u[(local_i - 1) * (BLOCK_Y + 2) + local_j] = d_u[(i - 1) * ny + j];
        sh_v[(local_i - 1) * (BLOCK_Y + 2) + local_j] = d_v[(i - 1) * (ny + 1) + j];
    }
    if (threadIdx.x == BLOCK_X - 1 && i < nx - 1) {
        sh_h[(local_i + 1) * (BLOCK_Y + 2) + local_j] = d_h[(i + 1) * (ny + 1) + j];
        sh_u[(local_i + 1) * (BLOCK_Y + 2) + local_j] = d_u[(i + 1) * ny + j];
        sh_v[(local_i + 1) * (BLOCK_Y + 2) + local_j] = d_v[(i + 1) * (ny + 1) + j];
    }
    if (threadIdx.y == 0 && j > 0) {
        sh_h[local_i * (BLOCK_Y + 2) + local_j - 1] = d_h[i * (ny + 1) + j - 1];
        sh_u[local_i * (BLOCK_Y + 2) + local_j - 1] = d_u[i * ny + j - 1];
        sh_v[local_i * (BLOCK_Y + 2) + local_j - 1] = d_v[i * (ny + 1) + j - 1];
    }
    if (threadIdx.y == BLOCK_Y - 1 && j < ny - 1) {
        sh_h[local_i * (BLOCK_Y + 2) + local_j + 1] = d_h[i * (ny + 1) + j + 1];
        sh_u[local_i * (BLOCK_Y + 2) + local_j + 1] = d_u[i * ny + j + 1];
        sh_v[local_i * (BLOCK_Y + 2) + local_j + 1] = d_v[i * (ny + 1) + j + 1];
    }

    __syncthreads();

    // Phase 1: Compute derivatives and store them in global memory
    if (i < nx - 1 && j < ny - 1) {
        double du_dx = (sh_u[(local_i + 1) * (BLOCK_Y + 2) + local_j] - sh_u[local_i * (BLOCK_Y + 2) + local_j]) / dx;
        double dv_dy = (sh_v[local_i * (BLOCK_Y + 2) + local_j + 1] - sh_v[local_i * (BLOCK_Y + 2) + local_j]) / dy;
        d_dh[i * ny + j] = -H * (du_dx + dv_dy);
    }

    if (i < nx - 1 && j < ny) {
        double dh_dx = (sh_h[(local_i + 1) * (BLOCK_Y + 2) + local_j] - sh_h[local_i * (BLOCK_Y + 2) + local_j]) / dx;
        d_du[i * ny + j] = -g * dh_dx;
    }

    if (i < nx && j < ny - 1) {
        double dh_dy = (sh_h[local_i * (BLOCK_Y + 2) + local_j + 1] - sh_h[local_i * (BLOCK_Y + 2) + local_j]) / dy;
        d_dv[i * ny + j] = -g * dh_dy;
    }

    __syncthreads();  // Ensure all derivatives are computed before proceeding to updates

    // Phase 2: Update fields using the derivatives calculated
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

// __global__ void compute_derivatives(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                                     double H, double g, double dx, double dy, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
//     int local_i = threadIdx.x + 1;  // Local index within shared memory, with padding for ghost cells
//     int local_j = threadIdx.y + 1;

//     extern __shared__ double shared_mem[];
//     double *sh_h = shared_mem;
//     double *sh_u = shared_mem + (BLOCK_X + 2) * (BLOCK_Y + 2);
//     double *sh_v = sh_u + (BLOCK_X + 2) * (BLOCK_Y + 2);

//     // Load h, u, v into shared memory with padding for ghost cells
//     if (i < nx && j < ny) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j] = d_h[i * (ny + 1) + j];
//         sh_u[local_i * (BLOCK_Y + 2) + local_j] = d_u[i * ny + j];
//         sh_v[local_i * (BLOCK_Y + 2) + local_j] = d_v[i * (ny + 1) + j];
//     }

//     // Load ghost cells for h, u, and v
//     if (threadIdx.x == 0 && i > 0) {
//         sh_h[(local_i - 1) * (BLOCK_Y + 2) + local_j] = d_h[(i - 1) * (ny + 1) + j];
//         sh_u[(local_i - 1) * (BLOCK_Y + 2) + local_j] = d_u[(i - 1) * ny + j];
//         sh_v[(local_i - 1) * (BLOCK_Y + 2) + local_j] = d_v[(i - 1) * (ny + 1) + j];
//     }
//     if (threadIdx.x == BLOCK_X - 1 && i < nx - 1) {
//         sh_h[(local_i + 1) * (BLOCK_Y + 2) + local_j] = d_h[(i + 1) * (ny + 1) + j];
//         sh_u[(local_i + 1) * (BLOCK_Y + 2) + local_j] = d_u[(i + 1) * ny + j];
//         sh_v[(local_i + 1) * (BLOCK_Y + 2) + local_j] = d_v[(i + 1) * (ny + 1) + j];
//     }
//     if (threadIdx.y == 0 && j > 0) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j - 1] = d_h[i * (ny + 1) + j - 1];
//         sh_u[local_i * (BLOCK_Y + 2) + local_j - 1] = d_u[i * ny + j - 1];
//         sh_v[local_i * (BLOCK_Y + 2) + local_j - 1] = d_v[i * (ny + 1) + j - 1];
//     }
//     if (threadIdx.y == BLOCK_Y - 1 && j < ny - 1) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j + 1] = d_h[i * (ny + 1) + j + 1];
//         sh_u[local_i * (BLOCK_Y + 2) + local_j + 1] = d_u[i * ny + j + 1];
//         sh_v[local_i * (BLOCK_Y + 2) + local_j + 1] = d_v[i * (ny + 1) + j + 1];
//     }

//     __syncthreads();

//     // Apply boundary conditions for h, u, v at the domain edges
//     if (i == nx) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j] = d_h[(nx - 1) * (ny + 1) + j];
//     }
//     if (j == ny) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j] = d_h[i * (ny + 1) + (ny - 1)];
//     }
//     if (i == 0) {
//         sh_u[local_i * BLOCK_Y + local_j] = d_u[(nx - 1) * ny + j];
//     }
//     if (j == 0) {
//         sh_v[local_i * (BLOCK_Y + 2) + local_j] = d_v[i * (ny + 1) + (ny - 1)];
//     }

//     __syncthreads();

//     // Compute derivatives using shared memory
//     if (i < nx - 1 && j < ny - 1) {
//         double du_dx = (sh_u[(local_i + 1) * (BLOCK_Y + 2) + local_j] - sh_u[local_i * (BLOCK_Y + 2) + local_j]) / dx;
//         double dv_dy = (sh_v[local_i * (BLOCK_Y + 2) + local_j + 1] - sh_v[local_i * (BLOCK_Y + 2) + local_j]) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }

//     if (i < nx - 1 && j < ny) {
//         double dh_dx = (sh_h[(local_i + 1) * (BLOCK_Y + 2) + local_j] - sh_h[local_i * (BLOCK_Y + 2) + local_j]) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }

//     if (i < nx && j < ny - 1) {
//         double dh_dy = (sh_h[local_i * (BLOCK_Y + 2) + local_j + 1] - sh_h[local_i * (BLOCK_Y + 2) + local_j]) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
// }


// __global__ void update_fields(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                               double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
//                               double a1, double a2, double a3, double dt, int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;

//     if (i >= nx || j >= ny) return;

//     // Update h, u, and v as in the original kernel
//     d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
//     if (i + 1 < nx) {
//         d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
//     }
//     if (j + 1 < ny) {
//         d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
//     }
// }

void step() {
    dim3 blockDim(BLOCK_X, BLOCK_Y);
    dim3 gridDim((nx + BLOCK_X - 1) / BLOCK_X, (ny + BLOCK_Y - 1) / BLOCK_Y);

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

    int shared_mem_size = (BLOCK_X + 2) * (BLOCK_Y + 2) * 3 * sizeof(double); // Shared memory size for h, u, v

    // Launch kernel to compute derivatives
    compute_and_update<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv,
                                   d_dh1, d_du1, d_dv1,
                                   d_dh2, d_du2, d_dv2,
                                   H,  g,  dx,  dy,
                                    a1, a2,  a3,  dt,
                                    nx,  ny);
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





