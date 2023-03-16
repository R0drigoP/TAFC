#include "molecula.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <ctime>
#include "mating_func.h"
#include "global.h"

using namespace std;

//global variables
int N_moleculas = 70;
int N_atomos = 3;
int dim_caixa = 10;
double survival_rate = 0.4;
double mutation_prob = 0.02;
int max_iter = 5000;
double **positions;

int parents_nb = int(survival_rate * N_moleculas + 0.5);
int couples_nb = int(parents_nb/2);
int children_per_couple = ( N_moleculas - parents_nb) / couples_nb;

bool mating = 1;

//probability of each molecule to be a parent
int main(){

  TCanvas *c1 = new TCanvas();
  auto gr = new TGraph();
  
  //population of molecules
  vector<molecula*> pop;

  for(int i = 0; i < N_moleculas; ++i)
    pop.push_back(new molecula(N_atomos, dim_caixa, mutation_prob));
 
  // is_parent[i] == 0 : False: Molecule i is not a parent
  // is_parent[i] == 1 : True: Molecule i is not a parent
  bool *is_parent = new bool[N_moleculas];
  int *parent_order = new int[couples_nb*2];

  for(int iter = 0; iter < max_iter; iter++){

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
      
      /*int survivors = int(survival_rate*N_moleculas);
      int replaced_per_surv = int((N_moleculas - survivors)/N_moleculas);
      
      cout << survivors + replaced_per_surv << endl;
      for(int i = 0; i < survivors; i++){
	for(int j = 0; j < replaced_per_surv; j++){
	  pop[survivors + j + i*replaced_per_surv] -> Set_Pos(pop[i] -> Get_Pos());
        }
      }
    }*/
    
      for(int mol=survival_rate*N_moleculas; mol<N_moleculas; mol+=survival_rate*N_moleculas){              
	for(int alive=0; alive<survival_rate*N_moleculas; ++alive)                          
	  pop[mol+alive]->Set_Pos(pop[alive]->Get_Pos());
      }
    }
    
    if( iter == max_iter-1){
      positions = pop[0] -> Get_Pos();
      for(int i = 0; i < N_atomos; i++)
	cout << "Atomo " << i << " : " << positions[i][0] << ", " << positions[i][1] << ", " << positions[i][2] << endl;
         
    }

  }

  delete[] is_parent;
  delete[] parent_order;

  pop.clear();

  c1 -> cd();
  gr -> Draw("AP");
  c1 -> SaveAs("evolution.pdf");
  
  return 0;
}
