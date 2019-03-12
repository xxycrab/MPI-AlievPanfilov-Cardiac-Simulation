/* 
 * Utilities for the Aliev-Panfilov code
 * Scott B. Baden, UCSD
 * Nov 2, 2015
 */

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

    /* First: allocate the sub-matrix for the process calling this function*/
    int nprocs, myrank;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    double *E_copy, *E_prev_copy, *R_copy;
    int sub_m = n/cb.py, sub_n = n/cb.px;  // initial sub-matrix size
    int r_m = n%(cb.py), r_n = n%(cb.px); // rest of row or column if each sub-matrix is of sub_m-by-sub_n
    if(myrank/cb.px < r_m) sub_m++;// due to the requirement of evenly distribute matrix
    if(myrank%cb.px < r_n) sub_n++;// due to the requirement of evenly distribute matrix

    // Allocate memory for the sub matrices
    mb.m = sub_m + 2;  // 2 more rows for ghost rows
    mb.n = sub_n + 2;  // 2 more columns for ghost columns
    mb.E_prev = (double*)malloc(mb.m*mb.n*sizeof(double));
    mb.R = (double*)malloc(mb.m*mb.n*sizeof(double));
    mb.E = (double*)malloc(mb.m*mb.n*sizeof(double));

    // Skip the first row which is filled with 0
    E_copy = E + (n+2) + 1;
    E_prev_copy = E_prev + (n+2) + 1;
    R_copy = R + (n+2) + 1;
    
    // This if branch is specific for process 0 because it does not receive data from another process
    if(myrank == 0) {
        for(int ii = 0 ; ii < sub_m ; ii+=1) {
            for(int jj = 0 ; jj < sub_n ; jj+=1) {
                mb.E_prev[(ii+1)*mb.n+jj+1] = E_prev_copy[ii*(n+2)+jj];
                mb.E[(ii+1)*mb.n+jj+1] = E_copy[ii*(n+2)+jj];
                mb.R[(ii+1)*mb.n+jj+1] = R_copy[ii*(n+2)+jj];
            }
        }
    }
    /* End of first step */

    /* Second: distribute jobs to other processors
     * For understanding, we assume a processor is assigned a submatrix of n/py * n/px
     * Then the processor matrix is py*px */
    int rows = n/cb.py, cols = n/cb.px; // initial sub-matrix size
    int rest_row = n%(cb.py), rest_col = n%(cb.px); // rest of row sor columns in global matrix if each sub-matrix is of rows-by-cols
    int incr_px = cb.px, incr_py = cb.py; // serve as the count for processors on a column (px) and on a row (py)
    incr_px--; // because process 0 has been dealt, decrease
    if(rest_row > 0) rows++; // due to the requirement of evenly distribute matrix
    if(rest_col > 1) cols++; // due to the requirement of evenly distribute matrix

    // Deal with the submatrix very next to the one for process 0.
    if(incr_px != 0) {
        E_prev = E_prev_copy + cols;
        R = R_copy + cols;
        E = E_copy + cols;
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
        else return;
    }

    // If it's process 0, it's responsible for distribute jobs
    if(myrank == 0) {
        for( int rank = 1 ; rank < nprocs ; rank++) {
            if(rank != 1) {
                incr_px--; // count down on number of processor on a row
                // first go through the processor on this row
                if(incr_px != 0) {
                    E_prev = E_prev + cols;
                    R = R + cols;
                    E = E + cols;
                }
                else {
                    incr_px = cb.px; // finish for a row of processor, restore thr count;
                    incr_py--; // count down on column
                    // begin with the first sub-matrix on the next processor row
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

            // due to the requirement of evenly distribution, increment rows or cols if there're rest row or column
            rows = n/cb.py;
            cols = n/cb.px;
            if(rest_row > cb.py - incr_py) rows++;
            if(rest_col > cb.px - incr_px) cols++;

            // allocate buffers for sending
            double* buffer_E_prev = (double*)malloc(rows*cols*sizeof(double));
            double* buffer_R = (double*)malloc(rows*cols*sizeof(double));
            double* buffer_E = (double*)malloc(rows*cols*sizeof(double));
            // fill buffer with data
            for(int ii = 0 ; ii < rows ; ii+=1) {
                for(int jj = 0 ; jj < cols ; jj+=1) {
                    buffer_E_prev[ii*cols+jj] = E_prev[(ii)*(n+2)+jj];
                    buffer_R[ii*cols+jj] = R[(ii)*(n+2)+jj];
                    buffer_E[ii*cols+jj] = E[(ii)*(n+2)+jj];
                }
            }
            // Send
            MPI_Send(buffer_E_prev,rows*cols,MPI_DOUBLE,rank,0,MPI_COMM_WORLD);
            MPI_Send(buffer_R,rows*cols,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
            MPI_Send(buffer_E,rows*cols,MPI_DOUBLE,rank,2,MPI_COMM_WORLD);
        }
    }
    // If not process 0, do not need to distribute jobs but need to receive the sub-matrices.
    else {
        // allocate buffers for receiving
        double* buffer_E_prev_tmp = (double*)malloc(sub_m*sub_n*sizeof(double));
        double* buffer_R_tmp = (double*)malloc(sub_m*sub_n*sizeof(double));
        double* buffer_E_tmp = (double*)malloc(sub_m*sub_n*sizeof(double));
        // Receive
        MPI_Recv(buffer_E_prev_tmp,sub_m*sub_n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(buffer_R_tmp,sub_m*sub_n,MPI_DOUBLE,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(buffer_E_tmp,sub_m*sub_n,MPI_DOUBLE,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        // Update myblock for solve function
        for(int ii = 0 ; ii < sub_m ; ii+=1) {
            for(int jj = 0 ; jj < sub_n ; jj+=1) {
                mb.E_prev[(ii+1)*mb.n+jj+1] = buffer_E_prev_tmp[ii*sub_n+jj];
                mb.R[(ii+1)*mb.n+jj+1] = buffer_R_tmp[ii*sub_n+jj];
                mb.E[(ii+1)*mb.n+jj+1] = buffer_E_tmp[ii*sub_n+jj];
            }
        }
    }
    /* End of second step */

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

double *alloc1D(int m, int n) {
    int nx = n, ny = m;
    double *E;
    // Ensures that allocatdd memory is aligned on a 16 byte boundary
    assert(E = (double *) memalign(16, sizeof(double) * nx * ny));
    return (E);
}


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
