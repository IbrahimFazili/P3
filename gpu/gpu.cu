// 1D, synch, 32
// Max error: 2.220446049250313e-14
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



// 2D, combined kernels, haloing, shared memory, synch, 16
// Max error: 2.220446049250313e-14
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

// __global__ void compute_and_update(double *d_h, double *d_u, double *d_v,
//                                    double *d_dh, double *d_du, double *d_dv,
//                                    double *d_dh1, double *d_du1, double *d_dv1,
//                                    double *d_dh2, double *d_du2, double *d_dv2,
//                                    double H, double g, double dx, double dy,
//                                    double a1, double a2, double a3, double dt,
//                                    int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
//     int local_i = threadIdx.x + 1;  // Local index within shared memory with padding for ghost cells
//     int local_j = threadIdx.y + 1;

//     extern __shared__ double shared_mem[];
//     double *sh_h = shared_mem;
//     double *sh_u = sh_h + (BLOCK_X + 2) * (BLOCK_Y + 2);
//     double *sh_v = sh_u + (BLOCK_X + 2) * (BLOCK_Y + 2);

//     // Load h, u, v into shared memory with padding for ghost cells
//     if (i < nx && j < ny) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j] = d_h[i * (ny + 1) + j];
//         sh_u[local_i * (BLOCK_Y + 2) + local_j] = d_u[i * ny + j];
//         sh_v[local_i * (BLOCK_Y + 2) + local_j] = d_v[i * (ny + 1) + j];
//     }

//     // Load halo regions for boundary values
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

//     // Phase 1: Compute derivatives and store them in global memory
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

//     __syncthreads();  // not needed // Ensure all derivatives are computed before proceeding to updates

//     // Phase 2: Update fields using the derivatives calculated
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
//     compute_and_update<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv,
//                                    d_dh1, d_du1, d_dv1,
//                                    d_dh2, d_du2, d_dv2,
//                                    H,  g,  dx,  dy,
//                                     a1, a2,  a3,  dt,
//                                     nx,  ny);
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



// 2D, combined kernels, haloing, shared memory, async, 16
// Max error: 2.220446049250313e-14
// Timing gpu size 100
// Total time: 1.021611
// Timing gpu size 500
// Total time: 1.2113299999999998
// Timing gpu size 1000
// Total time: 2.1425889999999996
// Timing gpu size 2000
// Total time: 6.434950000000001
// Timing gpu size 3000
// Total time: 12.02549
// Timing gpu size 4000
// Total time: 22.970968
// Timing gpu size 5000
// Total time: 33.767667
// Timing gpu size 6000
// Total time: 51.709748000000005
// Timing gpu size 7000
// Total time: 1.231563
// Timing gpu size 8000
// Total time: 92.46259
// Timing gpu size 9000
// Total time: 1.408479
// Timing gpu size 10000
// Total time: 143.93925900000002
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

// cudaStream_t stream;

// __global__ void compute_and_update(double *d_h, double *d_u, double *d_v,
//                                    double *d_dh, double *d_du, double *d_dv,
//                                    double *d_dh1, double *d_du1, double *d_dv1,
//                                    double *d_dh2, double *d_du2, double *d_dv2,
//                                    double H, double g, double dx, double dy,
//                                    double a1, double a2, double a3, double dt,
//                                    int nx, int ny) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
//     int local_i = threadIdx.x + 1;  // Local index within shared memory with padding for ghost cells
//     int local_j = threadIdx.y + 1;

//     extern __shared__ double shared_mem[];
//     double *sh_h = shared_mem;
//     double *sh_u = sh_h + (BLOCK_X + 2) * (BLOCK_Y + 2);
//     double *sh_v = sh_u + (BLOCK_X + 2) * (BLOCK_Y + 2);

//     // Load h, u, v into shared memory with padding for ghost cells
//     if (i < nx && j < ny) {
//         sh_h[local_i * (BLOCK_Y + 2) + local_j] = d_h[i * (ny + 1) + j];
//         sh_u[local_i * (BLOCK_Y + 2) + local_j] = d_u[i * ny + j];
//         sh_v[local_i * (BLOCK_Y + 2) + local_j] = d_v[i * (ny + 1) + j];
//     }

//     // Load halo regions for boundary values
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

//     // Phase 1: Compute derivatives and store them in global memory
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

//     __syncthreads();  // not needed // Ensure all derivatives are computed before proceeding to updates

//     // Phase 2: Update fields using the derivatives calculated
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

//     // Launch the kernel
//     compute_and_update<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv,
//                                    d_dh1, d_du1, d_dv1,
//                                    d_dh2, d_du2, d_dv2,
//                                    H, g, dx, dy,
//                                    a1, a2, a3, dt,
//                                    nx, ny);

//     // Synchronize to ensure kernel execution is complete
//     cudaDeviceSynchronize();

//     // Swap pointers for time-stepping
//     double *tmp;
//     tmp = d_dh2; d_dh2 = d_dh1; d_dh1 = d_dh; d_dh = tmp;
//     tmp = d_du2; d_du2 = d_du1; d_du1 = d_du; d_du = tmp;
//     tmp = d_dv2; d_dv2 = d_dv1; d_dv1 = d_dv; d_dv = tmp;

//     t++;
// }

// void transfer(double *h_host) {
//     cudaMemcpyAsync(h_host, d_h, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyDeviceToHost, stream);
//     cudaStreamSynchronize(stream);  // Ensure the transfer completes before accessing h_host
// }

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
//     cudaStreamDestroy(stream);  // Destroy the CUDA stream
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

//     // Allocate device memory
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

//     // Copy initial values to device
//     cudaMemcpy(d_h, h0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_u, u0, nx * ny * sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_v, v0, (nx + 1) * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice);

//     // Create the CUDA stream for asynchronous operations
//     cudaStreamCreate(&stream);
// }



// FASTEST SO FAR
// 1D, combined kernels, synch, 16
// Max error: 2.220446049250313e-14
// Timing gpu size 100
// Total time: 1.0088730000000001
// Timing gpu size 500
// Total time: 1.012065
// Timing gpu size 1000
// Total time: 1.7263680000000001
// Timing gpu size 2000
// Total time: 4.441633
// Timing gpu size 3000
// Total time: 8.937672
// Timing gpu size 4000
// Total time: 15.150151
// Timing gpu size 5000
// Total time: 23.200784
// Timing gpu size 6000
// Total time: 33.025794999999995
// Timing gpu size 7000
// Total time: 44.67675
// Timing gpu size 8000
// Total time: 57.91664600000001
// Timing gpu size 9000
// Total time: 73.271238
// Timing gpu size 10000
// Total time: 90.138533
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

#define BLOCKSIZE 512  // Adjusted for better occupancy on Perlmutter

__global__ void compute_and_multistep(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
                                      double *d_dh1, double *d_du1, double *d_dv1, double *d_dh2, double *d_du2, double *d_dv2,
                                      double H, double g, double dx, double dy, double a1, double a2, double a3, double dt,
                                      int nx, int ny) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i = idx / ny; // Convert linear idx to 2D grid coordinates
    int j = idx % ny;

    if (i >= nx || j >= ny) return; // Boundary check

    // Directly load values from global memory
    double h = d_h[i * (ny + 1) + j];
    double u = d_u[i * ny + j];
    double v = d_v[i * (ny + 1) + j];

    // Calculate ghost cells for h
    if (j < ny && i == nx) d_h[i * (ny + 1) + j] = d_h[(nx - 1) * (ny + 1) + j];  // Last row
    if (i < nx && j == ny) d_h[i * (ny + 1) + j] = d_h[i * (ny + 1) + (ny - 1)];  // Last column

    // Apply boundary conditions for u and v
    if (j < ny && i == 0) d_u[i * ny + j] = d_u[(nx - 1) * ny + j];  // First row
    if (i < nx && j == 0) d_v[i * (ny + 1) + j] = d_v[i * (ny + 1) + (ny - 1)]; // First column

    // Compute dh: finite differences along x and y using direct global memory access
    if (i < nx - 1 && j < ny - 1) {
        double du_dx = (d_u[(i + 1) * ny + j] - u) / dx;
        double dv_dy = (d_v[i * (ny + 1) + j + 1] - v) / dy;
        d_dh[i * ny + j] = -H * (du_dx + dv_dy);
    }

    // Compute du: finite differences along x using direct global memory access
    if (i < nx - 1 && j < ny) {
        double dh_dx = (d_h[(i + 1) * (ny + 1) + j] - h) / dx;
        d_du[i * ny + j] = -g * dh_dx;
    }

    // Compute dv: finite differences along y using direct global memory access
    if (i < nx && j < ny - 1) {
        double dh_dy = (d_h[i * (ny + 1) + j + 1] - h) / dy;
        d_dv[i * ny + j] = -g * dh_dy;
    }

    __syncthreads();

    // Multistep update for h, u, v
    d_h[i * (ny + 1) + j] += (a1 * d_dh[i * ny + j] + a2 * d_dh1[i * ny + j] + a3 * d_dh2[i * ny + j]) * dt;
    if (i + 1 < nx) {
        d_u[(i + 1) * ny + j] += (a1 * d_du[i * ny + j] + a2 * d_du1[i * ny + j] + a3 * d_du2[i * ny + j]) * dt;
    }
    if (j + 1 < ny) {
        d_v[i * (ny + 1) + j + 1] += (a1 * d_dv[i * ny + j] + a2 * d_dv1[i * ny + j] + a3 * d_dv2[i * ny + j]) * dt;
    }
}

void step() {
    dim3 gridDim((nx * ny + BLOCKSIZE - 1) / BLOCKSIZE); // Launch enough blocks to cover all elements
    dim3 blockDim(BLOCKSIZE);

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

    int shared_mem_size = 3 * BLOCKSIZE * sizeof(double);
    compute_and_multistep<<<gridDim, blockDim, shared_mem_size>>>(d_h, d_u, d_v, d_dh, d_du, d_dv,
                                                                  d_dh1, d_du1, d_dv1, d_dh2, d_du2, d_dv2,
                                                                  H, g, dx, dy, a1, a2, a3, dt, nx, ny);
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

