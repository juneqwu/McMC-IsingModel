/*Metropolis Hastings step intra-thread */

void MH(int size, double state[size][size], double new_state[size][size], double beta, int J, int k, int flip){
  int i,j;
  double uniform, like = 0;
  proposal(size, new_state, flip);
  like = exp(-beta*(hamiltonian(size, new_state, J, k) - hamiltonian(size, state,J,k)));
  uniform = rand() % 1000/1000.0;
  /*printf("%f \n", uniform);*/
  if (like > uniform){
    for (i = 0; i< size; i++){
      for (j = 0; j < size; j++){
	state[i][j] = new_state[i][j];
      }
    }
    /*printf("Proposal accepted! intra-thread!\n");*/
  }
  else{
     for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            new_state[i][j] = state[i][j];
         }
      }
     /* printf("Proposal rejected! intra-thread!\n");*/
  }
}

/* Metropolis Hastings step inter-thread */
void swap(int size, double myself[size*size], double mypartner[size*size], int J, int k, double beta){

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
  like = exp(beta*(hamiltonian(size, mypartner_m, J, k) - hamiltonian(size, myself_m,J,k)));
  /*printf("the likelihood is %f. \n", like);*/
  uniform = rand() % 1000/1000.0;
  /* printf("%f \n", uniform);*/
  if (like > uniform){
    for (i = 0; i< size*size; i++){
      temp[i] = myself[i];
      myself[i] = mypartner[i];
      mypartner[i] = temp[i];
      }
    /* printf("Proposal accepted! inter-thread!\n");*/
  }
}
