// error 12 with shared memory, 40 sec at 10,000
// Timing gpu size 100
// Total time: 1.198713
// Timing gpu size 500
// Total time: 1.12666
// Timing gpu size 1000
// Total time: 1.287851
// Timing gpu size 2000
// Total time: 2.3694729999999997
// Timing gpu size 3000
// Total time: 4.2054290000000005
// Timing gpu size 4000
// Total time: 6.887900999999999
// Timing gpu size 5000
// Total time: 10.393035999999999
// Timing gpu size 6000
// Total time: 14.584203
// Timing gpu size 7000
// Total time: 19.365472
// Timing gpu size 8000
// Total time: 25.184748
// Timing gpu size 9000
// Total time: 31.608581
// Timing gpu size 10000
// Total time: 38.806453
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

// // Define TILE_SIZE based on Perlmutter's shared memory limit
// #define TILE_SIZE 32

// // Combined kernel with shared memory for efficient data access
// template <const uint BLOCKSIZE>
// __global__ void compute_step(double *d_h, double *d_u, double *d_v, double *d_dh, double *d_du, double *d_dv,
//                              double H, double g, double dx, double dy, int nx, int ny) {

//     // Define shared memory for each field, with padding for boundary cells
//     __shared__ double sh_h[TILE_SIZE + 1][TILE_SIZE + 1];
//     __shared__ double sh_u[TILE_SIZE + 1][TILE_SIZE];
//     __shared__ double sh_v[TILE_SIZE][TILE_SIZE + 1];

//     // Calculate 2D indices within the block using a 1D thread layout
//     const int local_i = threadIdx.x / BLOCKSIZE;   // Row index within tile
//     const int local_j = threadIdx.x % BLOCKSIZE;   // Column index within tile

//     // Calculate global indices
//     const int i = blockIdx.x * TILE_SIZE + local_i;
//     const int j = blockIdx.y * TILE_SIZE + local_j;
    
//     // Load data into shared memory, with boundary checks
//     if (i <= nx && j <= ny) {
//         sh_h[local_i][local_j] = d_h[i * (ny + 1) + j];
//     }
//     if (i < nx && j < ny) {
//         sh_u[local_i][local_j] = d_u[i * ny + j];
//         sh_v[local_i][local_j] = d_v[i * (ny + 1) + j];
//     }

//     // Load extra data for ghost cells and boundaries
//     if (local_i == TILE_SIZE - 1 && i + 1 <= nx) {
//         sh_h[local_i + 1][local_j] = d_h[(i + 1) * (ny + 1) + j];
//     }
//     if (local_j == TILE_SIZE - 1 && j + 1 <= ny) {
//         sh_h[local_i][local_j + 1] = d_h[i * (ny + 1) + (j + 1)];
//     }

//     __syncthreads();

//     // Compute ghost cells for h
//     if (local_j < TILE_SIZE && local_i == TILE_SIZE - 1 && i == nx) {
//         sh_h[local_i + 1][local_j] = sh_h[0][local_j];
//     }
//     if (local_i < TILE_SIZE && local_j == TILE_SIZE - 1 && j == ny) {
//         sh_h[local_i][local_j + 1] = sh_h[local_i][0];
//     }

//     // Apply boundary conditions for u and v
//     if (j < ny && i == 0) {
//         sh_u[local_i][local_j] = sh_u[TILE_SIZE - 1][local_j];  // Horizontal boundary for u
//     }
//     if (i < nx && j == 0) {
//         sh_v[local_i][local_j] = sh_v[local_i][TILE_SIZE - 1];  // Vertical boundary for v
//     }

//     __syncthreads();

//     // Compute dh using shared memory
//     if (i < nx - 1 && j < ny - 1) {
//         double du_dx = (sh_u[local_i + 1][local_j] - sh_u[local_i][local_j]) / dx;
//         double dv_dy = (sh_v[local_i][local_j + 1] - sh_v[local_i][local_j]) / dy;
//         d_dh[i * ny + j] = -H * (du_dx + dv_dy);
//     }

//     // Compute du using shared memory
//     if (i < nx - 1 && j < ny) {
//         double dh_dx = (sh_h[local_i + 1][local_j] - sh_h[local_i][local_j]) / dx;
//         d_du[i * ny + j] = -g * dh_dx;
//     }

//     // Compute dv using shared memory
//     if (i < nx && j < ny - 1) {
//         double dh_dy = (sh_h[local_i][local_j + 1] - sh_h[local_i][local_j]) / dy;
//         d_dv[i * ny + j] = -g * dh_dy;
//     }
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
//     // dim3 gridDim(ceil(nx / TILE_SIZE), ceil(ny / TILE_SIZE));
//     // dim3 blockDim(TILE_SIZE, TILE_SIZE);
//     dim3 gridDim((nx + TILE_SIZE - 1) / TILE_SIZE, (ny + TILE_SIZE - 1) / TILE_SIZE);
//     dim3 blockDim(TILE_SIZE, TILE_SIZE);

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
//     dim3 blockDim(32, 32);

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

