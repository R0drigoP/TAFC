#include "molecula.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <ctime>
#include <fstream>
#include "mating_func.h"
#include "global.h"

using namespace std;

//global variables
int N_moleculas = 30;
int N_atomos = 13;
int dim_caixa = 10;
double survival_rate = 0.2;
double mutation_prob = 0.05;
int max_iter = 30000;

int parents_nb = int(survival_rate * N_moleculas + 0.5);
int couples_nb = int(parents_nb/2);
int children_per_couple = ( N_moleculas - parents_nb) / couples_nb;

bool mating = 1;

//probability of each molecule to be a parent
int main(){

  double **positions;

  if(mating == 0)
    survival_rate = 0.3;

  if(mating == 1)
    survival_rate = 0.45;

  TCanvas *c1 = new TCanvas();
  auto gr = new TGraph();

  ofstream text_file("best_molecule.bs");
  ofstream movie_file("best_molecule.mv");
  
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
    /*
    for(int i = 0; i < N_moleculas; i++){
      cout<<"Molecula"<<i<<endl;
      positions = pop[i] -> Get_Pos();
        cout << "Atomo " << 0 << " : " <<positions[0][0] << ", " << positions[0][1] << ", " << positions[0][2] << endl;
        cout << "Atomo " << 1 << " : " <<positions[1][0] << ", " << positions[1][1] << ", " << positions[1][2] << endl;
        cout << "Atomo " << 2 << " : " <<positions[2][0] << ", " << positions[2][1] << ", " << positions[2][2] << endl;

      }*/


    for(int i = 0; i < N_moleculas; ++i)
      pop[i] -> Potencial();

    //cout<<"got potential"<<endl;

    sort(pop.begin(), pop.end(), molecula::LessPot);

    //cout<<"sorted"<<endl;

    if(iter == 1){
      double** best_pos = pop[0]->Get_Pos();
      for(int i=0; i<N_atomos; ++i){
        text_file << "atom C " << flush;
        for(int j=0; j<3; ++j)
          text_file << best_pos[i][j] << " " << flush;
        text_file << endl;
      }
    }
    else if(iter%1000==0){
      movie_file << "frame" << endl;
      double** best_pos = pop[0]->Get_Pos();
      for(int i=0; i<N_atomos; ++i){
        for(int j=0; j<3; ++j)
          movie_file << best_pos[i][j] << " ";
      }
      movie_file << endl << endl;
    }

    gr -> AddPoint( iter, pop[0] -> Get_Pot());
    
    //cout<<"added point"<<endl;
    
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
      
      for(int mol = survival_rate*N_moleculas; mol < N_moleculas; mol += survival_rate*N_moleculas){              
        int alive = 0;
        while(alive < survival_rate*N_moleculas && (mol+alive)<N_moleculas){                         
          pop[mol+alive] -> Set_Pos(pop[alive] -> Get_Pos());
          ++alive;
        }
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
  
  for(int i = 0; i < 3; ++i) 
    delete[] positions[i];
  
  delete[] positions;
  
  pop.clear();
  
  text_file << endl << "spec C 0.1 Red" << endl
	    << endl << "bonds C C 0.5 1.5 0.01 0.0"
	    << endl << "bonds C H 0.4 1.0 0.01 1.0"
	    << endl << "scale 100"
	    << endl << "inc 3" << flush;
  
  c1 -> cd();
  gr -> Draw("AP");
  c1 -> SaveAs("evolution.pdf");
  
  text_file.close();
  movie_file.close();
  
  return 0;
}
