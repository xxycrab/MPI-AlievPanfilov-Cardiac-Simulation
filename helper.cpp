/* 
 * Utilities for the Aliev-Panfilov code
 * Scott B. Baden, UCSD
 * Nov 2, 2015
 */
//---
#include <iostream>
#include <assert.h>
// Needed for memalign
#include <malloc.h>
#include "cblock.h"

using namespace std;

typedef struct _my_block {

    int m, n;
    double *E;
    double *E_prev;
    double *R;

} my_block;

#ifdef _MPI_
#include <mpi.h>
#endif

extern control_block cb;

void printMat(const char mesg[], double *E, int m, int n);

my_block mb;

//
// Initialization
//
// We set the right half-plane of E_prev to 1.0, the left half plane to 0
// We set the botthom half-plane of R to 1.0, the top half plane to 0
// These coordinates are in world (global) coordinate and must
// be mapped to appropriate local indices when parallelizing the code
//
void init(double *E, double *E_prev, double *R, int m, int n) {
    int i;

    for (i = 0; i < (m + 2) * (n + 2); i++)
        E_prev[i] = R[i] = 0;

    for (i = (n + 2); i < (m + 1) * (n + 2); i++) {
        int colIndex = i % (n + 2);        // gives the base index (first row's) of the current index

        // Need to compute (n+1)/2 rather than n/2 to work with odd numbers
        if (colIndex == 0 || colIndex == (n + 1) || colIndex < ((n + 1) / 2 + 1))
            continue;

        E_prev[i] = 1.0;
    }

    for (i = 0; i < (m + 2) * (n + 2); i++) {
        int rowIndex = i / (n + 2);        // gives the current row number in 2D array representation
        int colIndex = i % (n + 2);        // gives the base index (first row's) of the current index

        // Need to compute (m+1)/2 rather than m/2 to work with odd numbers
        if (colIndex == 0 || colIndex == (n + 1) || rowIndex < ((m + 1) / 2 + 1))
            continue;

        R[i] = 1.0;
    }
    // We only print the meshes if they are small enough
#if 1
    printMat("E_prev", E_prev, m, n);
    printMat("R", R, m, n);
#endif

/************************************** Our code here ******************************************/
#ifdef _MPI_
    int nprocs, myrank;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    double *E_copy, *E_prev_copy, *R_copy;

    int sub_m = n/cb.py, sub_n = n/cb.px;  // matrix block size
    int r_m = n%(cb.py), r_n = n%(cb.px);
    if(myrank/cb.px < r_m) sub_m++;//?
    if(myrank%cb.px < r_n) sub_n++;//?

    // Allocate memory for a single matrix block
    mb.m = sub_m + 2;  // 2 more rows for ghost rows
    mb.n = sub_n + 2;  // 2 more columns for ghost columns
    mb.E_prev = (double*)malloc(mb.m*mb.n*sizeof(double));
    mb.R = (double*)malloc(mb.m*mb.n*sizeof(double));
    mb.E = (double*)malloc(mb.m*mb.n*sizeof(double));

    // Skip the first row which is filled with 0
    E_copy = E + (n+2) + 1;
    E_prev_copy = E_prev + (n+2) + 1;
    R_copy = R + (n+2) + 1;

    int rows = n/cb.py, cols = n/cb.px;
    int incr_row = n%(cb.py), incr_col = n%(cb.px);
    int incr_px = cb.px, incr_py = cb.py;
    incr_px--;
    if(incr_row > cb.py - incr_py) rows++;
    if(incr_col > cb.px - incr_px) cols++;

    // Copy elements to single matrix block
    if(myrank == 0) {
        for(int ii = 0 ; ii < sub_m ; ii+=1) {
            for(int jj = 0 ; jj < sub_n ; jj+=1) {
                mb.E_prev[(ii+1)*mb.n+jj+1] = E_prev_copy[ii*(n+2)+jj];
                mb.E[(ii+1)*mb.n+jj+1] = E_copy[ii*(n+2)+jj];
                mb.R[(ii+1)*mb.n+jj+1] = R_copy[ii*(n+2)+jj];
            }
        }
    }

    // 
    if(incr_px != 0) {
        E_prev = E_prev_copy + cols;
        R = R_copy + cols;
        E = E_copy + cols;
    }
    else {
        incr_px = cb.px;		//Added
        incr_py--;
        if(incr_py != 0) {
            E_copy = E_copy + rows*(n+2);
            E_prev_copy = E_prev_copy + rows*(n+2);
            R_copy = R_copy + rows*(n+2);
            E_prev = E_prev_copy;
            R = R_copy;
            E = E_copy;
        }
        else return;
    }

    if(myrank == 0) {
        for( int rank = 1 ; rank < nprocs ; rank++) {
            if(rank != 1) {
                incr_px--;
                if(incr_px != 0) {
                    E_prev = E_prev + cols;
                    R = R + cols;
                    E = E + cols;
                }
                else {
                    incr_px = cb.px;
                    incr_py--;
                    if(incr_py != 0) {
                        E_copy = E_copy + rows*(n+2);
                        E_prev_copy = E_prev_copy + rows*(n+2);
                        R_copy = R_copy + rows*(n+2);
                        E_prev = E_prev_copy;
                        R = R_copy;
                        E = E_copy;
                    }
                    else break;
                }
            }

            rows = n/cb.py;
            cols = n/cb.px;
            if(incr_row > cb.py - incr_py) rows++;
            if(incr_col > cb.px - incr_px) cols++;

            double* buffer_E_prev = (double*)malloc(rows*cols*sizeof(double));
            double* buffer_R = (double*)malloc(rows*cols*sizeof(double));
            double* buffer_E = (double*)malloc(rows*cols*sizeof(double));
            for(int ii = 0 ; ii < rows ; ii+=1) {
                for(int jj = 0 ; jj < cols ; jj+=1) {
                    buffer_E_prev[ii*cols+jj] = E_prev[(ii)*(n+2)+jj];
                    buffer_R[ii*cols+jj] = R[(ii)*(n+2)+jj];
                    buffer_E[ii*cols+jj] = E[(ii)*(n+2)+jj];
                }
            }
            MPI_Send(buffer_E_prev,rows*cols,MPI_DOUBLE,rank,0,MPI_COMM_WORLD);
            MPI_Send(buffer_R,rows*cols,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
            MPI_Send(buffer_E,rows*cols,MPI_DOUBLE,rank,2,MPI_COMM_WORLD);
        }
    }
    else {
        double* buffer_E_prev_tmp = (double*)malloc(sub_m*sub_n*sizeof(double));
        double* buffer_R_tmp = (double*)malloc(sub_m*sub_n*sizeof(double));
        double* buffer_E_tmp = (double*)malloc(sub_m*sub_n*sizeof(double));
        MPI_Recv(buffer_E_prev_tmp,sub_m*sub_n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(buffer_R_tmp,sub_m*sub_n,MPI_DOUBLE,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(buffer_E_tmp,sub_m*sub_n,MPI_DOUBLE,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        for(int ii = 0 ; ii < sub_m ; ii+=1) {
            for(int jj = 0 ; jj < sub_n ; jj+=1) {
                mb.E_prev[(ii+1)*mb.n+jj+1] = buffer_E_prev_tmp[ii*sub_n+jj];
                mb.R[(ii+1)*mb.n+jj+1] = buffer_R_tmp[ii*sub_n+jj];
                mb.E[(ii+1)*mb.n+jj+1] = buffer_E_tmp[ii*sub_n+jj];
            }
        }
    }
#else
    mb.m = m + 2;
    mb.n = n + 2;

    mb.E_prev = (double *) malloc(mb.m * mb.n * sizeof(double));
    mb.R = (double *) malloc(mb.m * mb.n * sizeof(double));
    mb.E = (double *) malloc(mb.m * mb.n * sizeof(double));

    mb.E_prev = E_prev;
    mb.E = E;
    mb.R = R;
#endif
/******************************************* End ***********************************************/
}

// No need to edit
double *alloc1D(int m, int n) {
    int nx = n, ny = m;
    double *E;
    // Ensures that allocatdd memory is aligned on a 16 byte boundary
    assert(E = (double *) memalign(16, sizeof(double) * nx * ny));
    return (E);
}


// No need to edit
void printMat(const char mesg[], double *E, int m, int n) {
    int i;
#if 0
    if (m>8)
        return;
#else
    if (m > 34)
        return;
#endif
    printf("%s\n", mesg);
    for (i = 0; i < (m + 2) * (n + 2); i++) {
        int rowIndex = i / (n + 2);
        int colIndex = i % (n + 2);
        if ((colIndex > 0) && (colIndex < n + 1))
            if ((rowIndex > 0) && (rowIndex < m + 1))
                printf("%6.3f ", E[i]);
        if (colIndex == n + 1)
            printf("\n");
    }
}
