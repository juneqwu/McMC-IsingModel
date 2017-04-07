/* computes the total magnetization of the system */

double magnetization(int size, double sigma[size][size]){
  int i,j;
  double magnetization = 0;

  for (i = 1; i < size-1; i++){
    for (j = 1; j < size-1; j++){
      magnetization += sigma[i][j];
    }
  }

  return magnetization;
}

/* computes the hamiltonian of the system */

double hamiltonian(int size, double sigma[size][size], double J, double k){
  int i,j;
  double hamiltonian = 0;
  for (i = 1; i < size-1; i++){
    for (j = 1; j < size-1; j++){
      hamiltonian += sigma[i][j] * (sigma[i-1][j-1]+ sigma[i-1][j] + sigma[i-1][j+1] +sigma[i+1][j-1]+ sigma[i+1][j+1]+ sigma[i+1][j] + sigma[i][j-1] + sigma[i][j+1]);
    }
  }

  double mag = magnetization(size, sigma);
  hamiltonian = -J*hamiltonian - k*mag*mag;
  return hamiltonian;
}
