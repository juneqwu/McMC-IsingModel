/* 2D Ising Model using Metropolis Hastings and Parallel Tempering.
   In parallel tempering algorithm, external magnetic field is used to
   predict critical temperature instead of using different levels of
   temperature in Boltzman distribution.

   Try it! The ising gets sweeter as the iterations go by!

   Author: June Wu */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "state.h"
#include "hamiltonian.h"
#include "MH.h"
#include <time.h>
#include "clcg4.h"


int accepts=0; /* Global variable to keep track of number of accepted swaps*/

int main(int argc, char * argv[]) {
  int nthreads, count;
  InitDefault(); /* Initialize random number generator with time*/


    /* Set the total number of threads from user input; The number
       threads corresponds to the number of parallel instance. */

    if (argc == 2) {
      nthreads = atol(argv[1]);
    } else {
      nthreads = 8;
    }
    
  printf("Welcome to the 2D Ising Model! \n");
  printf("The total number of parallel simulations is %d. \n", nthreads);
  
  for (count = 0; count <= 100; count++) {
    double T  = 0.5 + 3.5 * count / 100.0; /* Temperature range*/
    int L = 20; /* The grid size of the ising model. Total number of sites is L*L. */
    int size = L + 2; /* including the ghost points on the boundary */
    double J = 1, k = 0.000001, beta = 1.0 / (2.13*T), mag = 0; /* constants in the ising model */
    double state[size][size], new_state[size][size], ks[nthreads];
    int i, j, iter, partner, burn_in =1000, max_iter =1000000, m = 2500, flip = 1;

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
        MH(size, state, new_state, beta, J, ks[my_threadNum], flip);
      }

      #pragma omp barrier

      /* partner varialbe switches the left/right neighbour */
      partner = 1;

      /* Iterations with Swaps */
      for (iter = 1; iter <= max_iter; iter++) {

        /*if (my_threadNum == 0){
          printf("iteration %d, %f \n", iter, magnetization(size, state));
          } */

        /* Metropolis-Hasting step intra-thread */
 
        MH(size, state, new_state, beta, J, ks[my_threadNum], flip);
	mag += abs(magnetization(size, state));

        if (iter % m == 0){
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
	    swap(L, nthreads, my_threadNum, myself, mypartner, ks, J, ks[my_threadNum], beta, partner);
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
        }
      } /* end of iteration of swap */

        if (my_threadNum == 0 ) {
	   printf("Temperature = %f, Magnetization = %f\n", T, mag/max_iter);
	}
    } /* end of parallel routine */
            /* Write data to file*/
  } /* end of looping through temperature */
  printf("The total number of acceptance is %d.\n",accepts);
  return 0;
}
