/* 
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 */
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include "time.h"
#include "apf.h"
#include "cblock.h"
#include "Plotting.h"
#include <emmintrin.h>

#ifdef _MPI_
#include <mpi.h>
#endif
using namespace std;

typedef struct _my_block {
    int m, n;
    double *E;
    double *E_prev;
    double *R;
} my_block;


void repNorms(double l2norm, double mx, double dt, int m, int n, int niter, int stats_freq);

void stats(double *E, int m, int n, double *_mx, double *sumSq);

void printMat2(const char mesg[], double *E, int m, int n);

extern control_block cb;
extern my_block mb;

#ifdef SSE_VEC
// If you intend to vectorize using SSE instructions, you must
// disable the compiler's auto-vectorizer
__attribute__((optimize("no-tree-vectorize")))
#endif

// The L2 norm of an array is computed by taking sum of the squares
// of each element, normalizing by dividing by the number of points
// and then taking the sequare root of the result
//
double L2Norm(double sumSq) {
    double l2norm = sumSq / (double) ((cb.m) * (cb.n));
    l2norm = sqrt(l2norm);
    return l2norm;
}

void
solve(double **_E, double **_E_prev, double *R, double alpha, double dt, Plotter *plotter, double &L2, double &Linf) {

    // Simulated time is different from the integer timestep number
    double t = 0.0;
    //double *E = *_E, *E_prev = *_E_prev;
    double *E = mb.E, *E_prev = mb.E_prev;
    //double *R_tmp = R;
    double *R_tmp = mb.R;
    R = mb.R;
    //double *E_tmp = *_E;
    double *E_tmp = mb.E;
    //double *E_prev_tmp = *_E_prev;
    double *E_prev_tmp = mb.E_prev;
    double mx, sumSq, fsumSq, fLinf;
    int niter;
    //int m = cb.m, n=cb.n;
    int m = mb.m - 2, n = mb.n - 2;
    int innerBlockRowStartIndex = (n + 2) + 1;
    int innerBlockRowEndIndex = (((m + 2) * (n + 2) - 1) - (n)) - (n + 2);

    int rank = 0, np = 1;
#ifdef _MPI_
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
    int x1 = rank % cb.px, y1 = rank / cb.px;
    // We continue to sweep over the mesh until the simulation has reached
    // the desired number of iterations
    for (niter = 0; niter < cb.niters; niter++) {
        if (cb.debug && (niter == 0)) {
            stats(E_prev, m, n, &mx, &sumSq);
            double l2norm = L2Norm(sumSq);
            repNorms(l2norm, mx, dt, m, n, -1, cb.stats_freq);
            if (cb.plot_freq)
                plotter->updatePlot(E, -1, m + 1, n + 1);
        }

        /*
         * Copy data from boundary of the computational box to the
         * padding region, set up for differencing computational box's boundary
         *
         * These are physical boundary conditions, and are not to be confused
         * with ghost cells that we would use in an MPI implementation
         *
         * The reason why we copy boundary conditions is to avoid
         * computing single sided differences at the boundaries
         * which increase the running time of solve()
         *
         */

        // 4 FOR LOOPS set up the padding needed for the boundary conditions
        int i, j;

        // Fills in the TOP Ghost Cells
        if (y1 == 0) {
            for (i = 0; i < (n + 2); i++) {
                E_prev[i] = E_prev[i + (n + 2) * 2];
            }
        }
        // Fills in the RIGHT Ghost Cells
        if (x1 == cb.px - 1) {
            for (i = (n + 1); i < (m + 2) * (n + 2); i += (n + 2)) {
                E_prev[i] = E_prev[i - 2];
            }
        }
        // Fills in the LEFT Ghost Cells
        if (x1 == 0) {
            for (i = 0; i < (m + 2) * (n + 2); i += (n + 2)) {
                E_prev[i] = E_prev[i + 2];
            }
        }
        // Fills in the BOTTOM Ghost Cells
        if (y1 == cb.py - 1) {
            for (i = ((m + 2) * (n + 2) - (n + 2)); i < (m + 2) * (n + 2); i++) {
                E_prev[i] = E_prev[i - (n + 2) * 2];
            }
        }

//////////////////////////////////////////////////////////////////////////////

#define FUSED 0

/************************************* MPI part ****************************************/
#ifdef _MPI_
        double* buffer_right_sen = (double*)malloc(m*sizeof(double));
        double* buffer_left_sen = (double*)malloc(m*sizeof(double));
        double* buffer_right_rec = (double*)malloc(m*sizeof(double));
        double* buffer_left_rec = (double*)malloc(m*sizeof(double));

        if ((cb.px>1) || (cb.py>1)){
            // Put left column of the matrix in send left
            if (x1!=0){
                int j=0;
                for(i=n+3; i<(m+1)*(n+2); i+=n+2){
                    buffer_left_sen[j] =E_prev[i];
                    j++;
                }
           }
           //Put right column of matrix in send right
            if (x1!=cb.px-1){
                int j=0;
                for(i=n+2+n; i<(m+1)*(n+2); i+=n+2){
                    buffer_right_sen[j] =E_prev[i];
                    j++;
                }
            }

            //MPI sending and receiving communicator
            MPI_Request sendleft, sendright, sendtop, sendbot;
            MPI_Request recleft, recright, rectop, recbot;

            if(x1!=0 && !cb.noComm) {
                //left
                MPI_Isend(buffer_left_sen, m, MPI_DOUBLE, rank - 1, 0,MPI_COMM_WORLD , &(sendleft));
                MPI_Irecv(buffer_left_rec, m, MPI_DOUBLE, rank - 1, 0,MPI_COMM_WORLD , &(recleft));
            }
            if(x1!=cb.px-1 && !cb.noComm) {
               //right
               MPI_Isend(buffer_right_sen, m, MPI_DOUBLE, rank + 1, 0,MPI_COMM_WORLD , &(sendright));
               MPI_Irecv(buffer_right_rec, m, MPI_DOUBLE, rank + 1, 0,MPI_COMM_WORLD , &(recright));
            }
            if(y1!=0 && !cb.noComm) {
                // top
                MPI_Isend(E_prev+n+3, n, MPI_DOUBLE, rank - cb.px, 0,MPI_COMM_WORLD , &(sendtop));
                MPI_Irecv(E_prev+1, n, MPI_DOUBLE, rank - cb.px, 0,MPI_COMM_WORLD , &(rectop));
            }
            if(y1!=cb.py-1 && !cb.noComm) {
                // bottom
                MPI_Isend(E_prev+ m*(n+2)+1, n, MPI_DOUBLE, rank + cb.px, 0,MPI_COMM_WORLD , &(sendbot));
                MPI_Irecv(E_prev+ (m+1)*(n+2)+1, n, MPI_DOUBLE, rank + cb.px, 0,MPI_COMM_WORLD , &(recbot));
            }

            // Synchronization
            if(x1!=0 && !cb.noComm) {
                MPI_Wait(&(recleft),MPI_STATUS_IGNORE);
                MPI_Wait(&(sendleft),MPI_STATUS_IGNORE);
            }
            if(x1!=cb.px-1 && !cb.noComm) {
                MPI_Wait(&(recright),MPI_STATUS_IGNORE);
                MPI_Wait(&(sendright),MPI_STATUS_IGNORE);
            }
            if(y1!=0 && !cb.noComm) {
                MPI_Wait(&(rectop),MPI_STATUS_IGNORE);
                MPI_Wait(&(sendtop),MPI_STATUS_IGNORE);
            }

            if(y1!=cb.py-1 && !cb.noComm) {
                MPI_Wait(&(recbot),MPI_STATUS_IGNORE);
                MPI_Wait(&(sendbot),MPI_STATUS_IGNORE);
            }

            //Put left column of matrix in send left
            if (x1!=0){
                int j=0;
                for(i=n+2; i<(m+1)*(n+2); i+=n+2){
                    E_prev[i]=buffer_left_rec[j];
                    j++;
                }
            }
            //Put right column of matrix in send right
            if (x1!=cb.px-1){
                int j=0;
                for(i=n+3+n; i<(m+1)*(n+2); i+=n+2){
                    E_prev[i]=buffer_right_rec[j];
                    j++;
                }
            }
        }
#endif
/****************************************** End *********************************************/

/********************************* Computation part, to realize vectorization******************************************/
#ifdef FUSED
        // Solve for the excitation, a PDE
        for (j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j += (n + 2)) {
            E_tmp = E + j;
            E_prev_tmp = E_prev + j;
            R_tmp = R + j;
            for (i = 0; i < n; i++) {
                E_tmp[i] = E_prev_tmp[i] + alpha * (E_prev_tmp[i + 1] + E_prev_tmp[i - 1] - 4 * E_prev_tmp[i] +
                                                    E_prev_tmp[i + (n + 2)] + E_prev_tmp[i - (n + 2)]);
                E_tmp[i] += -dt *
                            (kk * E_prev_tmp[i] * (E_prev_tmp[i] - a) * (E_prev_tmp[i] - 1) + E_prev_tmp[i] * R_tmp[i]);
                R_tmp[i] += dt * (epsilon + M1 * R_tmp[i] / (E_prev_tmp[i] + M2)) *
                            (-R_tmp[i] - kk * E_prev_tmp[i] * (E_prev_tmp[i] - b - 1));
            }
        }
#else
        // Solve for the excitation, a PDE
//        for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
//            E_tmp = E + j;
//                E_prev_tmp = E_prev + j;
//                for(i = 0; i < n; i++) {
//                    E_tmp[i] = E_prev_tmp[i]+alpha*(E_prev_tmp[i+1]+E_prev_tmp[i-1]-4*E_prev_tmp[i]+E_prev_tmp[i+(n+2)]+E_prev_tmp[i-(n+2)]);
//                }
//        }
//
//        /*
//         * Solve the ODE, advancing excitation and recovery variables
//         *     to the next timtestep
//         */
//
//        for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
//            E_tmp = E + j;
//            R_tmp = R + j;
//        E_prev_tmp = E_prev + j;
//            for(i = 0; i < n; i++) {
//          E_tmp[i] += -dt*(kk*E_prev_tmp[i]*(E_prev_tmp[i]-a)*(E_prev_tmp[i]-1)+E_prev_tmp[i]*R_tmp[i]);
//          R_tmp[i] += dt*(epsilon+M1* R_tmp[i]/( E_prev_tmp[i]+M2))*(-R_tmp[i]-kk*E_prev_tmp[i]*(E_prev_tmp[i]-b-1));
//            }
//        }
        // use local variables instead of global to increase register hit rate
        double e_tmp, r_tmp;
        for (j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j += (n + 2)) {
            E_tmp = E + j;
            E_prev_tmp = E_prev + j;
            R_tmp = R + j;

            double e_tmps[n];
            // improve cache hit rate by put PDE and ODE in the same outer loop 
            // so that reduce cache missing when loading E_prev_tmp and E_tmp
            for (i = 0; i < n; i++) {
                e_tmps[i] = E_prev_tmp[i] + alpha * (E_prev_tmp[i + 1] + E_prev_tmp[i - 1] - 4 * E_prev_tmp[i] +
                                                     E_prev_tmp[i + (n + 2)] + E_prev_tmp[i - (n + 2)]);
            }
            for (i = 0; i < n; i++) {
                e_tmp = e_tmps[i];
                r_tmp = R_tmp[i];
                e_tmp += -dt * (kk * e_tmp * (e_tmp - a) * (e_tmp - 1) + e_tmp * r_tmp);
                r_tmp += dt * (epsilon + M1 * r_tmp / (e_tmp + M2)) * (-r_tmp - kk * e_tmp * (e_tmp - b - 1));
                E_tmp[i] = e_tmp;
                R_tmp[i] = r_tmp;
            }
        }

#endif

        /////////////////////////////////////////////////////////////////////////////////

        if (cb.stats_freq) {
            if (!(niter % cb.stats_freq)) {
                stats(E, m, n, &mx, &sumSq);
                double l2norm = L2Norm(sumSq);
                repNorms(l2norm, mx, dt, m, n, niter, cb.stats_freq);
            }
        }

        if (cb.plot_freq) {
            if (!(niter % cb.plot_freq)) {
                plotter->updatePlot(E, niter, m, n);
            }
        }

        // Swap current and previous meshes
        double *tmp = E;
        E = E_prev;
        E_prev = tmp;

    } //end of 'niter' loop at the beginning

//    printMat2("Rank 0 Matrix E_prev", E_prev, m, n); // return the L2 and infinity norms via in-out parameters

    stats(E_prev, m, n, &Linf, &sumSq);
/********************************************************************************************************************/

#ifdef _MPI_
    if (!cb.noComm){
        MPI_Reduce(&sumSq,&fsumSq ,1,MPI_DOUBLE,MPI_SUM,0 ,MPI_COMM_WORLD);
        MPI_Reduce(&Linf,&fLinf ,1,MPI_DOUBLE,MPI_MAX,0 ,MPI_COMM_WORLD);
        Linf = fLinf;
    }
    else fsumSq =sumSq;
#else
    fsumSq = sumSq;
#endif

    L2 = L2Norm(fsumSq);


#ifdef _MPI_
    free(mb.R);
    free(mb.E);
    free(mb.E_prev);
#else
    // Swap pointers so we can re-use the arrays
    *_E = E;
    *_E_prev = E_prev;
#endif
}

// No need to edit
void printMat2(const char mesg[], double *E, int m, int n) {
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
