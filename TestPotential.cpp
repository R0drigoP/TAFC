#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "molecula.h"

using namespace std;

int main()
{
  int N_moleculas = 10;
  int N_atomos = 10;
  int dim_caixa = 10;
  
  //Inicializar população como vector de ponteiros de objetos
  vector<molecula*> pop;

  //Preencher o vector
  for(int i=0; i<N_moleculas; ++i)
    pop.push_back(new molecula(N_atomos, dim_caixa));

  
  for(int i=0; i<N_moleculas; ++i){
    cout << "Juna: " << pop[i]->Potencial() << endl;
    cout << "Rodrigo: " << pop[i]->OtherPotential() << endl;
  }
  
  return 0;
}
