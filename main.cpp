#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "molecula.h"

using namespace std;

int main()
{
  int N_moleculas = 3;
  int N_atomos = 10;
  int dim_caixa = 10;
  
  double f_value;
  
  /*molecula* mol = new molecula(N_atomos, dim_caixa);
  
  f_value = mol -> Potencial();

  cout << f_value << endl;
  */

  vector<molecula*> pop;
  
  for(int i=0; i<N_moleculas; ++i)
    pop.push_back(new molecula(N_atomos, dim_caixa));

  
  for(int i=0; i<N_moleculas; ++i)
    pop[i]->Potencial();

  for(int i=0; i<N_moleculas; ++i)
    cout << pop[i]->Get_Fvalue() << endl;
  
  sort(pop.begin(), pop.end(), greater<molecula*>());

  cout << "====== Sorting Done ======" << endl;

  for(int i=0; i<N_moleculas; ++i)
    cout << pop[i]->Get_Fvalue() << endl;
  
  
  /*
  gRandom = new TRandom3(0);
  double x = gRandom->Uniform(-1.,1.);
  cout << x << endl;
  */
  
  return 0;
}
