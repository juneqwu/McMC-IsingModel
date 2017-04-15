/* 2D Ising Model using Metropolis simulation and Parallel Tempering.
   In parallel tempering algorithm, external mean field is used to 
   predict critical temperature instead of using different levels of
   temperature in Boltzman distribution. 

   Try it! And the ising gets sweeter as the iterations go by.
   (around 1e+6 in the cold temperature; 5e+6 around the hot 
    temperature).PT is a better sampler for the cold and the 
   critical temperature.

   Author: June Wu */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "state.h"
#include "hamiltonian.h"
#include "MH.h"
#include <time.h>

int main(int argc, char * argv[]) {
  int nthreads, count;
  double T = 1; /* Temperature */
  
    /* Set the total number of threads from user input; The number
       threads corresponds to the number of parallel instance. */
  
    if (argc == 2) {
      nthreads = atol(argv[1]);
    } else {
      nthreads = 8;
    }

  for (count = 0; count < 1; count++) {
    /*T  = 1 + 4.0 * count / 50.0;*/
    int L = 20; /* The grid size of the ising model. Total number of sites is L*L. */
    int size = L + 2; /* including the ghost points on the boudnary */
    double J = 1, k = 0.001, beta = 1.0 / (2.3*T), mag = 0; /* constancs in the ising model */
    double state[size][size], new_state[size][size], ks[nthreads];
    int i, j, iter, partner, burn_in =0, max_iter =50000, m = 1, flip = 1;

    /* A shared vector to help inter-thread swap */
    double neighbours[L * L][nthreads];

    omp_set_num_threads(nthreads);

#pragma omp parallel default (none) private(i, j, partner, iter, state, new_state, mag) shared(T, size, L, m, J, k, beta, flip, burn_in, max_iter, neighbours, ks, nthreads)
    {
      int my_threadNum = omp_get_thread_num();
      int numThreads = omp_get_num_threads();

      ks[my_threadNum] = k*my_threadNum;


      /* Initialize state and proposal state */
      initialize(size, state);
      for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
          new_state[i][j] = state[i][j];
        }
      }

      /* Burn-in iterations */
      for (i = 0; i < burn_in; i++) {
        MH(size, state, new_state, beta, J, k, flip);
      }

      #pragma omp barrier


      /* partner varialbe switches the left/right neighbour */
      partner = 1;

      /* Iterations with Swaps */
      for (iter = 0; iter < max_iter; iter++) {

	/*mag += hamiltonian(size, state, J, ks[my_threadNum]); */
	/*mag += abs(magnetization(size, state)); */
 
        /*if (my_threadNum == 0){
          printf("iteration %d, %f \n", iter, magnetization(size, state));
          } */

        /* Metropolis-Hasting step intra-thread */
        for (i = 0; i < m; i++) {
          MH(size, state, new_state, beta, J, k, flip);
        }

        /* Put states from all threads into a shared matrix neighbours */
        for (i = 1; i < size - 1; i++) {
          for (j = 1; j < size - 1; j++) {
            neighbours[(i - 1) * (size - 2) + (j - 1)][my_threadNum] = state[i][j];
          }
        }

        #pragma omp barrier
 
        /* Inter-thread Swap */
        #pragma omp critical
        if (my_threadNum % 2 == 1 && my_threadNum + partner < numThreads && my_threadNum + partner >= 0) {
          double myself[L * L], mypartner[L * L];
          for (i = 0; i < L * L; i++) {
            myself[i] = neighbours[i][my_threadNum];
            mypartner[i] = neighbours[i][my_threadNum + partner];
          }
	  swap(L, nthreads, my_threadNum, myself, mypartner, ks,J, k, beta, partner);
          for (i = 0; i < L * L; i++) {
            neighbours[i][my_threadNum] = myself[i];
            neighbours[i][my_threadNum + partner] = mypartner[i];
          }
        }
        partner = -partner;

        #pragma omp barrier

        /* Flush out swaped state to private matrix state */
        for (i = 1; i < size - 1; i++) {
          for (j = 1; j < size - 1; j++) {
            state[i][j] = neighbours[(i - 1) * (size - 2) + (j - 1)][my_threadNum];
          }
        }
	
       if (my_threadNum == 0 ) {
	 printf("%f,",  hamiltonian(size,state, J, ks[my_threadNum]));
      }

       
      } /* end of iteration of swap */
    } /* end of parallel routine */
  } /* end of looping through temperature */
  return 0;
}
