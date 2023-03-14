#include "molecula.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <ctime>
#include "mating_func.h"
#include "global.h"

using namespace std;

//global variables
int N_moleculas = 10;
int N_atomos = 3;
int dim_caixa = 10;
double survival_rate = 0.45;
double mutation_prob = 0.05;
int max_iter = 4;

int parents_nb = int(survival_rate * N_moleculas + 0.5);
int couples_nb = int(parents_nb/2);
int children_per_couple = ( N_moleculas - parents_nb) / couples_nb;

bool mating = 1;

//probability of each molecule to be a parent
int main(){

  double **positions;

  if(mating == 0)
    survival_rate = 0.2;

  if(mating == 1)
    survival_rate = 0.45;

  TCanvas *c1 = new TCanvas();
  auto gr = new TGraph();
  
  //population of molecules
  vector<molecula*> pop;

  positions = new double*[N_atomos];
  
  for (int i = 0; i < N_atomos; i++) 
    positions[i] = new double[3];

  for(int i = 0; i < N_moleculas; ++i)
    pop.push_back(new molecula(N_atomos, dim_caixa, mutation_prob));
 
  // is_parent[i] == 0 : False: Molecule i is not a parent
  // is_parent[i] == 1 : True: Molecule i is not a parent
  bool *is_parent = new bool[N_moleculas];
  int *parent_order = new int[couples_nb*2];

  for(int iter = 0; iter < max_iter; iter++){

    cout << "ITER NR " << iter << endl;

    for(int i = 0; i < N_moleculas; ++i)
      pop[i] -> OtherPotential();

    sort(pop.begin(), pop.end(), molecula::LessPot);

    gr -> AddPoint( iter, pop[0] -> Get_Pot());

    if( mating == 1){
      parent_probability( pop, is_parent, parent_order);
      generate_children( pop, parent_order);
    }

    if( mating == 0){

      for(int mol = 0; mol < N_moleculas; mol++){
      	pop[mol] -> Mutate();
      	pop[mol] -> OtherPotential();
      }
      
      sort(pop.begin(), pop.end(), molecula::LessPot);

      for(int mol = survival_rate*N_moleculas; mol < N_moleculas; mol += survival_rate*N_moleculas)              
        for(int alive = 0; alive < survival_rate*N_moleculas; ++alive)                          
          pop[mol+alive] -> Set_Pos(pop[alive] -> Get_Pos());
  }

  if( iter == max_iter-1){
    positions = pop[0] -> Get_Pos();
   
    for(int i = 0; i < N_atomos; i++)
      cout << "Atomo " << i << " : " << positions[i][0] << ", " << positions[i][1] << ", " << positions[i][2] << endl;
  }

  }

  delete[] is_parent;
  delete[] parent_order;

  for(int i = 0; i < 3; ++i) 
    delete[] positions[i];
    
  delete[] positions;

  pop.clear();

  c1 -> cd();
  gr -> Draw("AP");
  c1 -> SaveAs("evolution.pdf");
  
  return 0;
}
