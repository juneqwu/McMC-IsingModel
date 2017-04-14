/* initializes state */

void initialize(int size, double state[size][size]){
  int i, j;

  /* initialize boudnary ghost points as 0 */
  for (i = 0; i< size; i++){
    state[i][0] = 0;
    state[i][size-1] = 0;
  }
  for (j = 0; j < size; j++){
    state[0][j] = 0;
    state[size-1][j]=0;
  }

  /*initialize other points in the lattice randomly */
  for (i=1; i < size-1; i++){
    for (j = 1; j < size-1; j++){
      if (rand() % 10000 / 10000.0 < 0.5){
	state[i][j] = 1 ;
      }
      else
	state[i][j] = -1;
    }
  }
}

/* propose a new state */
void proposal(int size, double new_state[size][size], int flip){
  int i;
  for (i = 0; i < flip; i++){
    int m = rand() % (size-2) + 1;
    int n = rand() % (size-2) + 1;
    /*printf("flip (%d, %d).\n",m,n);*/
    new_state[m][n] = new_state[m][n]* (-1);
  }
}

/* print out state */
void print_state(int size, double state[size][size]){
    printf("\n");
    for (int i = 0; i <size; i++){
    for (int j = 0; j<size; j++){
      printf(" %2d ", (int) state[i][j]);
    }
    printf("\n");
  }
}

