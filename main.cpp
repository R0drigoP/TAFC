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
  double survival_rate = 0.4;
  
  //Inicializar população como vector de ponteiros de objetos
  vector<molecula*> pop;

  //Preencher o vector
  for(int i=0; i<N_moleculas; ++i)
    pop.push_back(new molecula(N_atomos, dim_caixa));

  
  for(int i=0; i<N_moleculas; ++i)
    cout << pop[i]->Potencial() << endl;

  //Faz a ordenação
  sort(pop.begin(), pop.end(), molecula::LessPot);

  cout << "====== Sorting Done ======" << endl;

  for(int i=0; i<N_moleculas; ++i)
    cout << pop[i]->Potencial() << endl;

  //cout << "====== Check if Sorting is well made ======" << endl;

  //for(int i=0; i<N_moleculas-1; ++i)
  //  cout << pop[i]->Potencial() - pop[i+1]->Potencial() << endl;

  //4 parents -> 2 couples
  // 3 children for each
  int parents_nb = int(survival_rate * N_moleculas + 0.5);
  int couples_nb = parents_nb/2;
  int children_nb = N_moleculas - parents_nb;
  int children_per_couple = children_nb / couples_nb;
  for(int i = 0; i<couples_nb;i++){
    for (int j = 0; j < children_per_couple; j++){
      double gene_prop = 1./(j+2);
      pop[parents_nb+i*children_per_couple+j]->Mating(pop[couples_nb*i],pop[couples_nb*i+1], gene_prop);
      //cout << parents_nb+i*children_per_couple+j << endl;
      //cout<<couples_nb*i<<" "<<couples_nb*i+1<<endl;



    }
  }  

    cout<<"NEXT GEN"<<endl;
    for(int i=0; i<N_moleculas; ++i)
      cout << pop[i]->Potencial() << endl;


  return 0;
}
