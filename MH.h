#include <time.h>
#include "clcg4.h"

/*Metropolis Hastings step intra-thread */

extern int accepts;

int MH(int size, double state[size][size], double new_state[size][size], double beta, double J, double k, int flip){
  int i,j;
  double uniform, like = 0;
  proposal(size, new_state, flip);
  like = exp(-beta*(hamiltonian(size, new_state, J, k) - hamiltonian(size, state,J,k)));
  /*printf("\nstate: %f \n", hamiltonian(size, state,J,k));
    printf("new state: %f \n", hamiltonian(size, new_state, J, k));*/
  uniform = GenVal(omp_get_thread_num());
  /*printf("random number %f with likelihood %f \n", uniform, like);*/
  if (like > uniform){
    for (i = 0; i< size; i++){
      for (j = 0; j < size; j++){
	state[i][j] = new_state[i][j];
      }
    }
    return 1;
    /*printf("Proposal accepted! intra-thread!\n");
      print_state(size, state);*/
  }
  else{
    /*printf("Proposal rejected! intra-thread!\n");
      print_state(size,new_state);*/
     for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            new_state[i][j] = state[i][j];
         }
      }
     return 0;
  }
  return 0;
}

/* Metropolis Hastings step inter-thread */
void swap(int size, int nthreads, int my_threadNum, double myself[size*size], double mypartner[size*size], double ks[nthreads], double J, double k, double beta, int partner){

  /*cast those two vectors into matrix so that it's easy to compute the hamiltonian */
  double myself_m[size][size], mypartner_m[size][size], temp[size*size];
  int i,j;
  double like, uniform;
  /* initialize boudnary ghost points as 0 */
  for (i = 0; i< size; i++){
    myself_m[i][0] = 0;
    myself_m[i][size-1] = 0;
    mypartner_m[i][0] = 0;
    mypartner_m[i][size-1] = 0;
  }
  for (j = 0; j < size; j++){
    myself_m[0][j] = 0;
    myself_m[size-1][j]=0;
    mypartner_m[0][j] = 0;
    mypartner_m[size-1][j]=0;
  }

  for (i = 0; i < size; i++){
    for (j = 0; j < size; j++){
      myself_m[i][j] = myself[i*size+j];
      mypartner_m[i][j] = mypartner[i*size+j];
    }
  }

  double delta_myself, delta_mypartner;
  delta_myself = hamiltonian(size, mypartner_m, J, k) - hamiltonian(size, myself_m, J, k);
  delta_mypartner = hamiltonian(size, mypartner_m, J, ks[my_threadNum + partner]) - hamiltonian(size, myself_m, J, ks[my_threadNum + partner]);

  like = exp(beta*(delta_mypartner - delta_myself));


  uniform = GenVal(my_threadNum);

  if (like > uniform){
    accepts +=1;
    for (i = 0; i< size*size; i++){
      temp[i] = myself[i];
      myself[i] = mypartner[i];
      mypartner[i] = temp[i];
      }
  }
}
