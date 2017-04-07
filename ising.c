#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "state.h"
#include "hamiltonian.h"
#include "MH.h"


int main(int argc, char *argv[]){
  double T = 1;
  for (int count = 0; count < 100; count++){
    T = 1 + 2*count/100;
  int nthreads, L = 20, size = L+2; /* including the ghost points on the boudnary */
  double J = 1, k = 0.1, beta = 1.0/T, mag = 0;
  double state[size][size], new_state[size][size];
  int i, j, iter, partner, burn_in = 500, max_iter = 1000, m = 1, flip = 10;

  /* set the total number of threads from user input */
  if (argc == 2){
    nthreads = atol(argv[1]);
  }
  else{
    nthreads = 8;
  }

  /* A shared vector to help inter-thread swap */
  double neighbours[L*L][nthreads];

  omp_set_num_threads(nthreads);

#pragma omp parallel default(none) private(i, j, partner, iter, state, new_state, mag) shared(size, L, m, J, k, beta, flip, burn_in, max_iter, neighbours)
  {
  int my_threadNum = omp_get_thread_num();
  int numThreads = omp_get_num_threads();

  k = k*my_threadNum;

  /* Initialize state and proposal state */
  initialize(size, state, my_threadNum);
  for (i = 0; i < size; i++){
    for (j = 0; j < size; j++){
       new_state[i][j] = state[i][j];
    }
  }


  /* Start the output file */
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "output%02d.txt", my_threadNum);
    fd = fopen(filename,"w+");
    fprintf(fd, "Thread %d out of %d says hi:\n", my_threadNum, numThreads);
    fprintf(fd, "The initial configuration is:\n");
    for(i = 0; i < size; i++){
      for(j = 0; j< size; j++){
	fprintf(fd, "%2d ", (int) state[i][j]);
      }
      fprintf(fd,"\n");
    }


    /* Burn-in iterations */
  for (i = 0; i < burn_in; i++){
    MH(size, state, new_state, beta, J, k, flip);
  }

  #pragma omp barrier

   /*printing results to a file*/
    fprintf(fd, "After burn-in configurations:\n");
    for(i = 0; i < size; i++){
      for(j = 0; j< size; j++){
	fprintf(fd, "%2d ", (int) state[i][j]);
      }
      fprintf(fd,"\n");
    }

    partner = 1;


      /* Iterations with Swaps */
  for (iter = 0; iter < max_iter; iter++){

    mag += magnetization(size, state);

    /* Metropolis-Hasting step intra-thread */
    for (i = 0; i < m; i++){
      MH(size, state, new_state, beta, J, k, flip);


      /*printing results to a file*/
      fprintf(fd, "Iteration %d with step size %d:\n",iter,m);
      for(i = 0; i < size; i++){
         for(j = 0; j< size; j++){
            fprintf(fd, "%2d ", (int) state[i][j]);
         }
         fprintf(fd,"\n");
      }

    }


    /* Put states from all threads into a shared matrix neighbours */
    for (i = 1; i < size-1; i++){
      for (j = 1; j < size-1; j++){
	neighbours[(i-1)*(size-2)+(j-1)][my_threadNum] = state[i][j];
      }
    }

#pragma omp barrier
       /*printing results to a file*/
    fprintf(fd, "shared info iteration %d after %d steps:\n",iter,m);
    int xi = 0;
    for (xi = 0; xi < numThreads; xi++){
      for (i = 0; i < L; i++){
	for (j = 0; j < L; j++){
	  fprintf(fd, "%2d ", (int) neighbours[i*L+j][xi]);
        }
       }
        fprintf(fd, "\n");
    }


    /* Inter-thread Swap */
#pragma omp critical
    if (my_threadNum %2 == 1 && my_threadNum + partner < numThreads && my_threadNum + partner >=0){
      double myself[L*L], mypartner[L*L];
      for (i = 0; i < L*L; i++){
	myself[i] = neighbours[i][my_threadNum];
	mypartner[i] = neighbours[i][my_threadNum + partner];
      }
      swap(L, myself, mypartner, J, k, beta);
      for (i = 0; i < L*L; i++){
        neighbours[i][my_threadNum] = myself[i];
        neighbours[i][my_threadNum + partner] = mypartner[i];
      }
    }
    partner = -partner;


#pragma omp barrier

    /*printing results to a file*/
    if (partner == -1)
      fprintf(fd, "After the swap with the right neighbour!\n");
    else
      fprintf(fd, "After the swap with the left neighbour!\n");
    for (xi = 0; xi < numThreads; xi++){
      for (i = 0; i < L; i++){
	for (j = 0; j < L; j++){
	  fprintf(fd, "%2d ", (int) neighbours[i*L+j][xi]);
        }
       }
        fprintf(fd, "\n");
    }

     /* Flush out swaped state to private matrix state */
 for (i = 1 ; i < size -1; i++){
   for (j = 1 ; j < size -1; j++){
     state[i][j] = neighbours[(i-1)*(size-2)+(j-1)][my_threadNum];
   }
   }
  }
  if (my_threadNum == 0){
    printf("%f,", mag/(max_iter*m));
  }
  fclose(fd);
  }
  }
  return 0;
}

