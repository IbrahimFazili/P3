// ibrahim work
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>

// #include "../common/common.hpp"
// #include "../common/solver.hpp"
// #define idx(i, j) ((i) * (local_ny + 1) + (j))

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
//     v = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * global_ny;
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h + (local_ny + 3), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(u0, sendcounts, displs, MPI_DOUBLE, u + (local_ny + 3), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v + (local_ny + 3), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h + (local_ny + 3), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, u + (local_ny + 3), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v + (local_ny + 3), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;
//     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;

//     for (int j = 1; j <= local_ny; j++) {
//         send_buffer[j - 1] = h[idx(1, j)];
//     }
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 0, recv_buffer, local_ny, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h[idx(local_nx + 1, j)] = recv_buffer[j - 1];
//     }

//     for (int j = 1; j <= local_ny; j++) {
//         send_buffer[j - 1] = h[idx(local_nx, j)];
//     }
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 1, recv_buffer, local_ny, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h[idx(0, j)] = recv_buffer[j - 1];
//     }
// }

// void compute_derivatives()
// {
//     for (int i = 1; i <= local_nx; i++) {
//         for (int j = 1; j <= local_ny; j++) {
//             double dx_h = (h[idx(i + 1, j)] - h[idx(i - 1, j)]) / (2 * dx);
//             double dy_h = (h[idx(i, j + 1)] - h[idx(i, j - 1)]) / (2 * dy);
            
//             double dx_u = (u[idx(i + 1, j)] - u[idx(i - 1, j)]) / (2 * dx);
//             double dy_v = (v[idx(i, j + 1)] - v[idx(i, j - 1)]) / (2 * dy);

//             dh[idx(i - 1, j - 1)] = -H * (dx_u + dy_v);
//             du[idx(i - 1, j - 1)] = -g * dx_h;
//             dv[idx(i - 1, j - 1)] = -g * dy_h;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }

//     for (int i = 1; i <= local_nx; i++) {
//         for (int j = 1; j <= local_ny; j++) {
//             h[idx(i, j)] += (a1 * dh[idx(i - 1, j - 1)] + a2 * dh1[idx(i - 1, j - 1)] + a3 * dh2[idx(i - 1, j - 1)]) * dt;
//             u[idx(i, j)] += (a1 * du[idx(i - 1, j - 1)] + a2 * du1[idx(i - 1, j - 1)] + a3 * du2[idx(i - 1, j - 1)]) * dt;
//             v[idx(i, j)] += (a1 * dv[idx(i - 1, j - 1)] + a2 * dv1[idx(i - 1, j - 1)] + a3 * dv2[idx(i - 1, j - 1)]) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * global_ny;
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h + local_ny + 3, local_nx * local_ny, MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }




// w macros
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
//     v = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * global_ny;
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, &h(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(u0, sendcounts, displs, MPI_DOUBLE, &u(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, &v(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, &h(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, &u(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, &v(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;
//     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;

//     for (int j = 1; j <= local_ny; j++) {
//         send_buffer[j - 1] = h(1, j);
//     }
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 0, recv_buffer, local_ny, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(local_nx + 1, j) = recv_buffer[j - 1];
//     }

//     for (int j = 1; j <= local_ny; j++) {
//         send_buffer[j - 1] = h(local_nx, j);
//     }
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 1, recv_buffer, local_ny, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(0, j) = recv_buffer[j - 1];
//     }
// }

// void compute_derivatives()
// {
//     for (int i = 1; i <= local_nx; i++) {
//         for (int j = 1; j <= local_ny; j++) {
//             double dx_h = (h(i + 1, j) - h(i - 1, j)) / (2 * dx);
//             double dy_h = (h(i, j + 1) - h(i, j - 1)) / (2 * dy);
            
//             double dx_u = (u(i + 1, j) - u(i - 1, j)) / (2 * dx);
//             double dy_v = (v(i, j + 1) - v(i, j - 1)) / (2 * dy);

//             dh(i - 1, j - 1) = -H * (dx_u + dy_v);
//             du(i - 1, j - 1) = -g * dx_h;
//             dv(i - 1, j - 1) = -g * dy_h;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }

//     for (int i = 1; i <= local_nx; i++) {
//         for (int j = 1; j <= local_ny; j++) {
//             h(i, j) += (a1 * dh(i - 1, j - 1) + a2 * dh1(i - 1, j - 1) + a3 * dh2(i - 1, j - 1)) * dt;
//             u(i, j) += (a1 * du(i - 1, j - 1) + a2 * du1(i - 1, j - 1) + a3 * du2(i - 1, j - 1)) * dt;
//             v(i, j) += (a1 * dv(i - 1, j - 1) + a2 * dv1(i - 1, j - 1) + a3 * dv2(i - 1, j - 1)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * global_ny;
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(&h(1, 1), local_nx * local_ny, MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }


// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx) * (local_ny + 2), sizeof(double));
//     v = (double*)calloc((local_nx) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * global_ny;
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, &h(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(u0, sendcounts, displs, MPI_DOUBLE, &u(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, &v(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, &h(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, &u(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, &v(1, 1), local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int row = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
//     int halo = (rank == 0) ? MPI_PROC_NULL : rank - 1;

//     // for (int j = 1; j <= local_ny; j++) {
//     std::memcpy(send_buffer, h + local_ny)
//     send_buffer[j - 1] = h(1, j);
//     // }
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, halo, 0, recv_buffer, local_ny, MPI_DOUBLE, row, 0, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(local_nx + 1, j) = recv_buffer[j - 1];
//     }

//     for (int j = 1; j <= local_ny; j++) {
//         send_buffer[j - 1] = h(local_nx, j);
//     }
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, row, 1, recv_buffer, local_ny, MPI_DOUBLE, halo, 1, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(0, j) = recv_buffer[j - 1];
//     }
// }


// void compute_derivatives()
// {
//     for (int i = 1; i <= local_nx; i++) {
//         for (int j = 1; j <= local_ny; j++) {
//             double dx_h = (h(i + 1, j) - h(i - 1, j)) / (2 * dx);
//             double dy_h = (h(i, j + 1) - h(i, j - 1)) / (2 * dy);
            
//             double dx_u = (u(i + 1, j) - u(i - 1, j)) / (2 * dx);
//             double dy_v = (v(i, j + 1) - v(i, j - 1)) / (2 * dy);

//             dh(i - 1, j - 1) = -H * (dx_u + dy_v);
//             du(i - 1, j - 1) = -g * dx_h;
//             dv(i - 1, j - 1) = -g * dy_h;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }

//     for (int i = 1; i <= local_nx; i++) {
//         for (int j = 1; j <= local_ny; j++) {
//             h(i, j) += (a1 * dh(i - 1, j - 1) + a2 * dh1(i - 1, j - 1) + a3 * dh2(i - 1, j - 1)) * dt;
//             u(i, j) += (a1 * du(i - 1, j - 1) + a2 * du1(i - 1, j - 1) + a3 * du2(i - 1, j - 1)) * dt;
//             v(i, j) += (a1 * dv(i - 1, j - 1) + a2 * dv1(i - 1, j - 1) + a3 * dv2(i - 1, j - 1)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * global_ny;
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(&h(1, 1), local_nx * local_ny, MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }



// -15
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx) * (local_ny + 1), sizeof(double));
//     u = (double*)calloc((local_nx) * local_ny, sizeof(double));
//     v = (double*)calloc((local_nx) * (local_ny + 1), sizeof(double));
    
//     dh = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = nx_i * global_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
//     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;

//     // std::memcpy(send_buffer, h + local_nx * local_ny, local_ny * sizeof(double));

//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 0, recv_buffer, local_ny, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
    
//     // for (int j = 1; j <= local_ny; j++) {
//     //     h(local_nx + 1, j) = recv_buffer[j - 1];
//     // }

// }

// void compute_derivatives()
// {
//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             double dhdx = (h(i + 1, j) - h(i, j)) / dx;
//             double dhdy = (h(i, j + 1) - h(i, j)) / dy;

//             dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
//             du(i, j) = -g * dhdx;
//             dv(i, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }

//??
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx) * (local_ny + 1), sizeof(double));
//     u = (double*)calloc((local_nx) * local_ny, sizeof(double));
//     v = (double*)calloc((local_nx) * (local_ny + 1), sizeof(double));
    
//     dh = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = nx_i * global_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
//     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;

//     // std::memcpy(send_buffer, h + local_nx * local_ny, local_ny * sizeof(double));

//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 0, recv_buffer, local_ny, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
    
//     // for (int j = 1; j <= local_ny; j++) {
//     //     h(local_nx + 1, j) = recv_buffer[j - 1];
//     // }

// }

// void compute_derivatives()
// {
//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             double dhdx = (h(i + 1, j) - h(i, j)) / dx;
//             double dhdy = (h(i, j + 1) - h(i, j)) / dy;

//             dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
//             du(i, j) = -g * dhdx;
//             dv(i, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }


//kaushik
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     v = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     du = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     dv = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny+2, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny+2, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = nx_i * global_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
//     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;

//     std::memcpy(send_buffer, h + local_nx * local_ny, local_ny * sizeof(double));

//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 0, recv_buffer, local_ny, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
    
//     for (int j = 1; j <= local_ny; j++) {
//         h(local_nx + 1, j) = recv_buffer[j - 1];
//     }

// }

// void compute_derivatives()
// {
//     for (int i = 1; i < local_nx+1; i++)
//     {
//         for (int j = 1; j < local_ny+1; j++)
//         {
//             double dhdx = (h(i , j) - h(i-1, j)) / dx;
//             double dhdy = (h(i-1, j + 1) - h(i-1, j)) / dy;

//             dh(i-1, j) = -H * (du_dx(i-1, j) + dv_dy(i-1, j));
//             du(i-1, j) = -g * dhdx;
//             dv(i-1, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }


//max error 19.479626705617992
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     v = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     du = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     dv = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny+2, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny+2, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = nx_i * global_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
//     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;

//     // Send the last row to the right and receive from the left
//     std::memcpy(send_buffer, h + local_nx * (local_ny + 1), local_ny * sizeof(double));
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 0, recv_buffer, local_ny, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(0, j) = recv_buffer[j - 1];
//     }

//     // Send the first row to the left and receive from the right
//     std::memcpy(send_buffer, h + 1 * (local_ny + 1), local_ny * sizeof(double));
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 0, recv_buffer, local_ny, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(local_nx + 1, j) = recv_buffer[j - 1];
//     }
// }

// void compute_derivatives()
// {
//     for (int i = 1; i < local_nx+1; i++)
//     {
//         for (int j = 1; j < local_ny+1; j++)
//         {
//             double dhdx = (h(i , j) - h(i-1, j)) / dx;
//             double dhdy = (h(i - 1, j + 1) - h(i - 1, j)) / dy;

//             dh(i-1, j) = -H * (du_dx(i-1, j) + dv_dy(i - 1, j));
//             du(i-1, j) = -g * dhdx;
//             dv(i-1, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }


//also 19
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     v = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     du = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     dv = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny+2, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny+2, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = nx_i * global_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, (local_nx+2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD); //changed this
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, (local_nx+2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
//     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;

//     // Send the last row to the right and receive from the left
//     std::memcpy(send_buffer, h + local_nx * (local_ny + 1), local_ny * sizeof(double));
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 0, recv_buffer, local_ny, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(0, j) = recv_buffer[j - 1];
//     }

//     // Send the first row to the left and receive from the right
//     std::memcpy(send_buffer, (h + 1) * (local_ny + 1), local_ny * sizeof(double));
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 0, recv_buffer, local_ny, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(local_nx + 1, j) = recv_buffer[j - 1];
//     }
// }

// void compute_derivatives()
// {
//     for (int i = 1; i < local_nx+1; i++)
//     {
//         for (int j = 1; j < local_ny+1; j++)
//         {
//             double dhdx = (h(i , j) - h(i - 1, j)) / dx;
//             double dhdy = (h(i - 1, j + 1) - h(i - 1, j)) / dy;

//             dh(i-1, j) = -H * (du_dx(i-1, j) + dv_dy(i - 1, j));
//             du(i-1, j) = -g * dhdx;
//             dv(i-1, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }

//wip
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     v = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     du = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     dv = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny+2, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny+2, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = nx_i * global_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, (local_nx+2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD); //changed this
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, (local_nx+2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// // void exchange_ghost_cells()
// // {
// //     MPI_Status status;
// //     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
// //     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    
// //     // Exchange ghost cells along the right and left boundaries
// //     // Copy last row to send buffer for right neighbor
// //     std::memcpy(send_buffer, &h(local_nx, 0), local_ny * sizeof(double));
// //     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 1, recv_buffer, local_ny, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, &status);
// //     for (int j = 0; j < local_ny; j++) {
// //         h(0, j) = recv_buffer[j];
// //     }
    
// //     // Copy first row to send buffer for left neighbor
// //     std::memcpy(send_buffer, &h(1, 0), local_ny * sizeof(double));
// //     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 2, recv_buffer, local_ny, MPI_DOUBLE, right, 2, MPI_COMM_WORLD, &status);
// //     for (int j = 0; j < local_ny; j++) {
// //         h(local_nx + 1, j) = recv_buffer[j];
// //     }

// //     // Add top-bottom exchanges similarly if needed.
// // }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
//     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;
//     int up = (rank == 0) ? MPI_PROC_NULL : rank - 1;
//     int down = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
    
//     // Buffers for top-bottom exchanges (need size `local_nx + 2`)
//     double *send_buffer_tb = (double*)calloc(local_nx + 2, sizeof(double));
//     double *recv_buffer_tb = (double*)calloc(local_nx + 2, sizeof(double));

//     // Exchange ghost cells along the right and left boundaries
//     // Send the last column to the right and receive from the left
//     std::memcpy(send_buffer, &h(local_nx, 0), local_ny * sizeof(double));
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 1, recv_buffer, local_ny, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, &status);
//     for (int j = 0; j < local_ny; j++) {
//         h(0, j) = recv_buffer[j];
//     }
    
//     // Send the first column to the left and receive from the right
//     std::memcpy(send_buffer, &h(1, 0), local_ny * sizeof(double));
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 2, recv_buffer, local_ny, MPI_DOUBLE, right, 2, MPI_COMM_WORLD, &status);
//     for (int j = 0; j < local_ny; j++) {
//         h(local_nx + 1, j) = recv_buffer[j];
//     }

//     // Top-Bottom Exchanges
//     // Send the last row to the bottom and receive from the top
//     std::memcpy(send_buffer_tb, &h(0, local_ny), (local_nx + 2) * sizeof(double));
//     MPI_Sendrecv(send_buffer_tb, local_nx + 2, MPI_DOUBLE, down, 3, recv_buffer_tb, local_nx + 2, MPI_DOUBLE, up, 3, MPI_COMM_WORLD, &status);
//     for (int i = 0; i < local_nx + 2; i++) {
//         h(i, 0) = recv_buffer_tb[i];
//     }
    
//     // Send the first row to the top and receive from the bottom
//     std::memcpy(send_buffer_tb, &h(0, 1), (local_nx + 2) * sizeof(double));
//     MPI_Sendrecv(send_buffer_tb, local_nx + 2, MPI_DOUBLE, up, 4, recv_buffer_tb, local_nx + 2, MPI_DOUBLE, down, 4, MPI_COMM_WORLD, &status);
//     for (int i = 0; i < local_nx + 2; i++) {
//         h(i, local_ny + 1) = recv_buffer_tb[i];
//     }

//     // Free the top-bottom exchange buffers
//     free(send_buffer_tb);
//     free(recv_buffer_tb);
// }


// void compute_derivatives()
// {
//     for (int i = 1; i < local_nx+1; i++)
//     {
//         for (int j = 1; j < local_ny+1; j++)
//         {
//             double dhdx = (h(i , j) - h(i - 1, j)) / dx;
//             double dhdy = (h(i - 1, j + 1) - h(i - 1, j)) / dy;

//             dh(i-1, j) = -H * (du_dx(i-1, j) + dv_dy(i - 1, j));
//             du(i-1, j) = -g * dhdx;
//             dv(i-1, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }


// 19.479626705617992 w/ block haloing
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     v = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     du = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     dv = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny+2, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny+2, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = nx_i * global_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, (local_nx+2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD); //changed this
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, (local_nx+2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int right = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;
//     int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;
//     int up = (rank == 0) ? MPI_PROC_NULL : rank - 1;
//     int down = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;

//     // left-right exchanges
//     // send the last column to the right and receive from the left
//     std::memcpy(send_buffer, h + local_nx * (local_ny + 1), local_ny * sizeof(double));
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, right, 0, recv_buffer, local_ny, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(0, j) = recv_buffer[j - 1];
//     }

//     // send the first column to the left and receive from the right
//     std::memcpy(send_buffer, h + (1 * (local_ny + 1)), local_ny * sizeof(double));
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, left, 0, recv_buffer, local_ny, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
//     for (int j = 1; j <= local_ny; j++) {
//         h(local_nx + 1, j) = recv_buffer[j - 1];
//     }

//     // top-bottom exchanges
//     // send the last row to the bottom and receive from the top
//     std::memcpy(send_buffer, h + (local_nx * (local_ny + 1)), (local_ny + 1) * sizeof(double));  // copy last row
//     MPI_Sendrecv(send_buffer, local_ny + 1, MPI_DOUBLE, down, 1, recv_buffer, local_ny + 1, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &status);
//     for (int i = 1; i <= local_nx; i++) {
//         h(i, 0) = recv_buffer[i - 1];
//     }

//     // send the first row to the top and receive from the bottom
//     std::memcpy(send_buffer, h + (1 * (local_ny + 1)), (local_ny + 1) * sizeof(double));  // copy first row
//     MPI_Sendrecv(send_buffer, local_ny + 1, MPI_DOUBLE, up, 1, recv_buffer, local_ny + 1, MPI_DOUBLE, down, 1, MPI_COMM_WORLD, &status);
//     for (int i = 1; i <= local_nx; i++) {
//         h(i, local_ny + 1) = recv_buffer[i - 1];
//     }
// }


// void compute_derivatives()
// {
//     for (int i = 1; i < local_nx+1; i++)
//     {
//         for (int j = 1; j < local_ny+1; j++)
//         {
//             double dhdx = (h(i , j) - h(i - 1, j)) / dx;
//             double dhdy = (h(i - 1, j + 1) - h(i - 1, j)) / dy;

//             dh(i-1, j) = -H * (du_dx(i-1, j) + dv_dy(i - 1, j));
//             du(i-1, j) = -g * dhdx;
//             dv(i-1, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }




//20, row only haloing
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     v = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc((local_nx+2) * (local_ny + 2), sizeof(double));
//     du = (double*)calloc((local_nx+2) * (local_ny+1), sizeof(double));
//     dv = (double*)calloc((local_nx+1) * (local_ny + 2), sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny+2, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny+2, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = nx_i * global_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, (local_nx+2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD); //changed this
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, (local_nx+2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int up = (rank == 0) ? MPI_PROC_NULL : rank - 1;
//     int down = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;

//     // top-bottom exchanges
//     // send the last row to the bottom and receive from the top
//     std::memcpy(send_buffer, h + (local_nx * (local_ny + 1)), (local_ny + 1) * sizeof(double));  // copy last row
//     MPI_Sendrecv(send_buffer, local_ny + 1, MPI_DOUBLE, down, 1, recv_buffer, local_ny + 1, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &status);
//     for (int i = 1; i <= local_nx; i++) {
//         h(i, 0) = recv_buffer[i - 1];
//     }

//     // send the first row to the top and receive from the bottom
//     std::memcpy(send_buffer, h + (1 * (local_ny + 1)), (local_ny + 1) * sizeof(double));  // copy first row
//     MPI_Sendrecv(send_buffer, local_ny + 1, MPI_DOUBLE, up, 1, recv_buffer, local_ny + 1, MPI_DOUBLE, down, 1, MPI_COMM_WORLD, &status);
//     for (int i = 1; i <= local_nx; i++) {
//         h(i, local_ny + 1) = recv_buffer[i - 1];
//     }
// }


// void compute_derivatives()
// {
//     for (int i = 1; i < local_nx+1; i++)
//     {
//         for (int j = 1; j < local_ny+1; j++)
//         {
//             double dhdx = (h(i , j) - h(i - 1, j)) / dx;
//             double dhdy = (h(i - 1, j + 1) - h(i - 1, j)) / dy;

//             dh(i-1, j) = -H * (du_dx(i-1, j) + dv_dy(i - 1, j));
//             du(i-1, j) = -g * dhdx;
//             dv(i-1, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }

// ed disussion cpy
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
//     u = (double*)calloc((local_nx + 2) * (local_ny + 1), sizeof(double));
//     v = (double*)calloc((local_nx + 1) * (local_ny + 2), sizeof(double));
    
//     dh = (double*)calloc((local_nx + 2) * (local_ny + 2), sizeof(double));
//     du = (double*)calloc((local_nx + 2) * (local_ny + 1), sizeof(double));
//     dv = (double*)calloc((local_nx + 1) * (local_ny + 2), sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny + 2, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny + 2, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = nx_i * global_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, (local_nx + 2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, (local_nx+2) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 2) * (local_ny+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, (local_nx +1) * (local_ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int up = (rank == 0) ? num_procs - 1 : rank - 1;
//     int down = (rank == num_procs - 1) ? 0 : rank + 1;

//     // top-bottom exchanges
//     // send the last row to the bottom and receive from the top
//     std::memcpy(send_buffer, h + (local_nx * (local_ny + 1)), (local_ny + 1) * sizeof(double));  // copy last row
//     MPI_Sendrecv(send_buffer, local_ny + 1, MPI_DOUBLE, down, 1, recv_buffer, local_ny + 1, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &status);
//     for (int i = 1; i <= local_nx; i++) {
//         h(i, 0) = recv_buffer[i - 1];
//     }

//     // send the first row to the top and receive from the bottom
//     std::memcpy(send_buffer, h + (1 * (local_ny + 1)), (local_ny + 1) * sizeof(double));  // copy first row
//     MPI_Sendrecv(send_buffer, local_ny + 1, MPI_DOUBLE, up, 1, recv_buffer, local_ny + 1, MPI_DOUBLE, down, 1, MPI_COMM_WORLD, &status);
//     for (int i = 1; i <= local_nx; i++) {
//         h(i, local_ny + 1) = recv_buffer[i - 1];
//     }
// }


// void compute_derivatives()
// {
//     for (int i = 1; i < local_nx+1; i++)
//     {
//         for (int j = 1; j < local_ny+1; j++)
//         {
//             double dhdx = (h(i , j) - h(i - 1, j)) / dx;
//             double dhdy = (h(i - 1, j + 1) - h(i - 1, j)) / dy;

//             dh(i - 1, j) = -H * (du_dx(i-1, j) + dv_dy(i - 1, j));
//             du(i - 1, j) = -g * dhdx;
//             dv(i - 1, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * (global_ny + 1);
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, local_nx * (local_ny + 1), MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }



//ours
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx + 1) * local_ny, sizeof(double));
//     u = (double*)calloc((local_nx + 1) * local_ny, sizeof(double));
//     v = (double*)calloc(local_nx * (local_ny + 1), sizeof(double));
    
//     dh = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         sendcounts[i] = local_nx * local_ny;  // Calculate the number of elements for each process
//         displs[i] = offset;               // Assign starting index for each process's section
//         offset += sendcounts[i];          // Update offset for the next process
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         u_sendcounts[i] = (local_nx + 1) * local_ny;  // Calculate the number of elements for each process
//         u_displs[i] = u_offset;               // Assign starting index for each process's section
//         u_offset += u_sendcounts[i];          // Update offset for the next process
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, (local_nx + 1) * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 1) * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, (local_nx + 1) * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 1) * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, local_nx * (local_ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
    
//     delete[] u_sendcounts;
//     delete[] u_displs;
// }

// void exchange_ghost_cells()
// {
//     MPI_Status status;
//     int up = (rank == 0) ? num_procs - 1 : rank - 1;
//     int down = (rank == num_procs - 1) ? 0 : rank + 1; 

//     // top-bottom exchanges
//     // send the last row to the bottom and receive from the top
//     std::memcpy(send_buffer, h + (local_nx - 1) * local_ny, local_ny * sizeof(double));  // copy last row
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, down, 1, recv_buffer, local_ny, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &status);
//     for (int i = 0; i < local_ny; i++) {
//         h(0, i) = recv_buffer[i];
//     }
// }

// void compute_derivatives()
// {
//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             double dhdx = (h(i + 1, j) - h(i, j)) / dx;
//             double dhdy = (h(i, j + 1) - h(i, j)) / dy;

//             dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
//             du(i, j) = -g * dhdx;
//             dv(i, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 0; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = local_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = (nx_i + 1) * global_ny;
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h, (local_nx + 1) * local_ny, MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }




//wip
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     // Allocate h and u with one extra row for top halo
//     h = (double*)calloc((local_nx + 1) * local_ny, sizeof(double));
//     u = (double*)calloc((local_nx + 1) * local_ny, sizeof(double));
    
//     // Allocate v with one extra row for top halo (no extra columns needed)
//     v = (double*)calloc(local_nx * (local_ny + 1), sizeof(double));
    
//     // Auxiliary arrays without halos
//     dh = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     // Communication buffers
//     send_buffer = (double*)calloc(local_ny, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     // Configuration for MPI_Scatterv
//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         sendcounts[i] = nx_i * local_ny;
//         displs[i] = offset;
//         offset += sendcounts[i];
//     }

//     int *u_sendcounts = new int[num_procs];
//     int *u_displs = new int[num_procs];
//     int u_offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         u_sendcounts[i] = (nx_i + 1) * local_ny;
//         u_displs[i] = u_offset;
//         u_offset += u_sendcounts[i];
//     }

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h, (local_nx + 1) * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
//         MPI_Scatterv(u0, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 1) * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h, (local_nx + 1) * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
//         MPI_Scatterv(nullptr, u_sendcounts, u_displs, MPI_DOUBLE, u, (local_nx + 1) * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;
//     delete[] u_sendcounts;
//     delete[] u_displs;

//     printf("Rank %d - local_nx: %d, local_ny: %d, global_nx: %d, global_ny: %d\n", rank, local_nx, local_ny, global_nx, global_ny);
//     printf("Rank %d - h size: %lu, u size: %lu, v size: %lu\n", rank, (local_nx + 1) * local_ny, (local_nx + 1) * local_ny, (local_nx + 1) * local_ny);
//     fflush(stdout);
// }


// void exchange_ghost_cells() {
//     MPI_Status status;
//     int up = (rank == 0) ? num_procs - 1 : rank - 1;
//     int down = (rank == num_procs - 1) ? 0 : rank + 1;

//     // Copy last row into send_buffer for sending
//     std::memcpy(send_buffer, h + (local_nx - 1) * local_ny, local_ny * sizeof(double));
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, down, 1, recv_buffer, local_ny, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &status);

    // // Place received data into the top halo row of h
    // for (int i = 0; i < local_ny; i++) {
    //     h[i] = recv_buffer[i];  // Correcting index to match the halo row location
    // }
// }

// void compute_derivatives()
// {
//     // #pragma omp parallel for collapse(2)
//     for (int i = 1; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             double dhdx = (h(i + 1, j) - h(i, j)) / dx;
//             double dhdy = (h(i, j + 1) - h(i, j)) / dy;

//             dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
//             du(i, j) = -g * dhdx;
//             dv(i, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 1; i < local_nx; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         // Calculate the local rows each process contributes, excluding halo rows
//         int nx_i = global_nx / num_procs + (i < global_nx % num_procs ? 1 : 0);
//         recvcounts[i] = nx_i * global_ny;  // Only core rows (no halos)
//         displs[i] = offset;
//         offset += recvcounts[i];
//     }

//     // Gather only the main data (excluding halos) from each process
//     MPI_Gatherv(h, local_nx * local_ny, MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     delete[] recvcounts;
//     delete[] displs;
// }


// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }



// max error 12
// #include <mpi.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <cstring>

// #define ny local_ny
// #include "../common/common.hpp"
// #include "../common/solver.hpp"

// int rank, num_procs;
// int local_nx, local_ny;
// int global_nx, global_ny;

// double *h, *u, *v;
// double *dh, *du, *dv;
// double *dh1, *du1, *dv1;
// double *dh2, *du2, *dv2;
// double *send_buffer, *recv_buffer;

// double H, g, dx, dy, dt;
// int t = 0;

// void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
// {
//     rank = rank_;
//     num_procs = num_procs_;
//     global_nx = nx_;
//     global_ny = ny_;
    
//     local_nx = global_nx / num_procs;
//     if (rank < global_nx % num_procs) {
//         local_nx++;
//     }
//     local_ny = global_ny;

//     h = (double*)calloc((local_nx + 2) * local_ny, sizeof(double));
//     u = (double*)calloc((local_nx + 1) * local_ny, sizeof(double));
//     v = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
//     dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
//     dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

//     send_buffer = (double*)calloc(local_ny, sizeof(double));
//     recv_buffer = (double*)calloc(local_ny, sizeof(double));

//     H = H_;
//     g = g_;
//     dx = length_ / global_nx;
//     dy = width_ / global_ny;
//     dt = dt_;

//     int *sendcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int local_nx_i = global_nx / num_procs;
//         if (i < global_nx % num_procs) {
//             local_nx_i++;  // Processes with an extra row
//         }
//         sendcounts[i] = local_nx_i * global_ny;  // Total elements for each process (local_nx * local_ny)
//         displs[i] = offset;                      // Starting position in the global array
//         offset += sendcounts[i];                 // Update offset for the next process
//     }

//     printf("Rank %d - local_nx: %d, local_ny: %d, global_nx: %d, global_ny: %d\n", rank, local_nx, local_ny, global_nx, global_ny);
//     fflush(stdout);
//     printf("Rank %d - h size: %lu, u size: %lu, v size: %lu\n", rank, (local_nx + 2) * local_ny, (local_nx + 1) * local_ny, local_nx * (local_ny + 1));
//     fflush(stdout);

//     if (rank == 0) {
//         MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h + local_ny, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
//         MPI_Scatterv(u0, sendcounts, displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h + local_ny, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     }

//     delete[] sendcounts;
//     delete[] displs;

// }

// void exchange_ghost_cells() {
//     MPI_Status status;
//     int up = (rank == 0) ? num_procs - 1 : rank - 1;
//     int down = (rank == num_procs - 1) ? 0 : rank + 1; 

//     // Send the last data row to the bottom and receive the top ghost row
//     std::memcpy(send_buffer, h + (local_nx - 1) * local_ny, local_ny * sizeof(double));  // Copy last data row
//     MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, down, 1, recv_buffer, local_ny, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &status);
//     std::memcpy(h, recv_buffer, local_ny * sizeof(double));  // Place into the top ghost row
// }

// void compute_derivatives()
// {
//     for (int i = 1; i < local_nx - 1; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             double dhdx = (h(i + 1, j) - h(i, j)) / dx;
//             double dhdy = (h(i, j + 1) - h(i, j)) / dy;

//             dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
//             du(i, j) = -g * dhdx;
//             dv(i, j) = -g * dhdy;
//         }
//     }
// }

// void step()
// {
//     exchange_ghost_cells();
//     compute_derivatives();

//     double a1, a2, a3;
//     if (t == 0) { a1 = 1.0; a2 = a3 = 0.0; }
//     else if (t == 1) { a1 = 3.0 / 2.0; a2 = -1.0 / 2.0; a3 = 0.0; }
//     else { a1 = 23.0 / 12.0; a2 = -16.0 / 12.0; a3 = 5.0 / 12.0; }


//     for (int i = 1; i < local_nx - 1; i++)
//     {
//         for (int j = 0; j < local_ny; j++)
//         {
//             h(i, j) += (a1 * dh(i, j) + a2 * dh1(i, j) + a3 * dh2(i, j)) * dt;
//             u(i + 1, j) += (a1 * du(i, j) + a2 * du1(i, j) + a3 * du2(i, j)) * dt;
//             v(i, j + 1) += (a1 * dv(i, j) + a2 * dv1(i, j) + a3 * dv2(i, j)) * dt;
//         }
//     }

//     double *tmp;
//     tmp = dh2; dh2 = dh1; dh1 = dh; dh = tmp;
//     tmp = du2; du2 = du1; du1 = du; du = tmp;
//     tmp = dv2; dv2 = dv1; dv1 = dv; dv = tmp;

//     t++;
// }

// void transfer(double *h_recv)
// {
//     int *recvcounts = new int[num_procs];
//     int *displs = new int[num_procs];
//     int offset = 0;

//     for (int i = 0; i < num_procs; i++) {
//         int local_nx_i = global_nx / num_procs;
//         if (i < global_nx % num_procs) {
//             local_nx_i++;  // Extra row for processes with uneven division
//         }
//         recvcounts[i] = local_nx_i * global_ny;  // Only the data portion
//         displs[i] = offset;                      // Offset in the global array
//         offset += recvcounts[i];
//     }

//     MPI_Gatherv(h + local_ny, local_nx * local_ny, MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

// }

// void free_memory()
// {
//     free(h); free(u); free(v);
//     free(dh); free(du); free(dv);
//     free(dh1); free(du1); free(dv1);
//     free(dh2); free(du2); free(dv2);
//     free(send_buffer);
//     free(recv_buffer);
// }





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

    h = (double*)calloc((local_nx + 2) * local_ny, sizeof(double));
    u = (double*)calloc((local_nx + 1) * local_ny, sizeof(double));
    v = (double*)calloc(local_nx * local_ny, sizeof(double));
    
    dh = (double*)calloc(local_nx * local_ny, sizeof(double));
    du = (double*)calloc(local_nx * local_ny, sizeof(double));
    dv = (double*)calloc(local_nx * local_ny, sizeof(double));
    
    dh1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    du1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    dv1 = (double*)calloc(local_nx * local_ny, sizeof(double));
    
    dh2 = (double*)calloc(local_nx * local_ny, sizeof(double));
    du2 = (double*)calloc(local_nx * local_ny, sizeof(double));
    dv2 = (double*)calloc(local_nx * local_ny, sizeof(double));

    send_buffer = (double*)calloc(local_ny, sizeof(double));
    recv_buffer = (double*)calloc(local_ny, sizeof(double));

    H = H_;
    g = g_;
    dx = length_ / global_nx;
    dy = width_ / global_ny;
    dt = dt_;

    int *sendcounts = new int[num_procs];
    int *displs = new int[num_procs];
    int offset = 0;

    for (int i = 0; i < num_procs; i++) {
        int local_nx_i = global_nx / num_procs;
        if (i < global_nx % num_procs) {
            local_nx_i++;  // Processes with an extra row
        }
        sendcounts[i] = local_nx_i * global_ny;  // Total elements for each process (local_nx * local_ny)
        displs[i] = offset;                      // Starting position in the global array
        offset += sendcounts[i];                 // Update offset for the next process
    }

    if (rank == 0) {
        MPI_Scatterv(h0, sendcounts, displs, MPI_DOUBLE, h + local_ny, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        MPI_Scatterv(u0, sendcounts, displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(v0, sendcounts, displs, MPI_DOUBLE, v, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, h + local_ny, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, u, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, sendcounts, displs, MPI_DOUBLE, v, local_nx * local_ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    delete[] sendcounts;
    delete[] displs;

}

void exchange_ghost_cells() {
    MPI_Status status;
    int up = (rank == 0) ? num_procs - 1 : rank - 1;
    int down = (rank == num_procs - 1) ? 0 : rank + 1; 

    // Send the last data row to the bottom and receive the top ghost row
    std::memcpy(send_buffer, h + (local_nx - 1) * local_ny, local_ny * sizeof(double));  // Copy last data row
    MPI_Sendrecv(send_buffer, local_ny, MPI_DOUBLE, down, 1, recv_buffer, local_ny, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &status);
    std::memcpy(h, recv_buffer, local_ny * sizeof(double));  // Place into the top ghost row
}

void compute_derivatives()
{
    for (int i = 1; i < local_nx - 1; i++)
    {
        for (int j = 0; j < local_ny; j++)
        {
            double dhdx = (h(i + 1, j) - h(i, j)) / dx;
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


    for (int i = 1; i < local_nx - 1; i++)
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
        int local_nx_i = global_nx / num_procs;
        if (i < global_nx % num_procs) {
            local_nx_i++;  // Extra row for processes with uneven division
        }
        recvcounts[i] = local_nx_i * global_ny;  // Only the data portion
        displs[i] = offset;                      // Offset in the global array
        offset += recvcounts[i];
    }

    MPI_Gatherv(h + local_ny, local_nx * local_ny, MPI_DOUBLE, h_recv, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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