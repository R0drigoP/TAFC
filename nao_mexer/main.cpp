#include "molecula.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <ctime>

using namespace std;

//global variables
int N_moleculas = 100;
int N_atomos = 10;
int dim_caixa = 3;
double survival_rate = 0.4;
double mutation_prob = 0.02;
int max_iter = 40000;

int parents_nb = int(survival_rate * N_moleculas + 0.5);
int couples_nb = int(parents_nb/2);
int children_per_couple = ( N_moleculas - parents_nb) / couples_nb;


//probability of each molecule to be a parent
void parent_probability(vector<molecula*> pop, bool *is_parent, int *parent_order ){
  
  double prob[N_moleculas], soma_prob_new[N_moleculas];
  double soma_prob = 0, soma_prob_aux = 0;

  // is_parent[i] == 0 : False: Molecule i is not a parent
  // is_parent[i] == 1 : True: Molecule i is not a parent

  int nr_couples_formed = 0, parent = 0;
  
  double r;
  srand((unsigned int)time(NULL));

  for(int i = 0; i < N_moleculas; i++){ 
    is_parent[i] = 0;
    soma_prob_new[i] = 0;
  }

  for(int i = 0; i < N_moleculas; i++){
    prob[i] =  exp( -(pop[i] -> Get_Pot()) );
    soma_prob += prob[i]; 
  }

  for(int i = 0; i < N_moleculas; i++){
    prob[i] = prob[i]/soma_prob;
    soma_prob_new[i] += soma_prob_aux; 
    soma_prob_aux += prob[i];
  }

    //select couples of parents (de acordo c/ pagina 76/77 Evolutionary Optimization Algorithms)
    //!! Acho que os pais deviam ser replaced pelas suas children (Bibliografia: Evolutionary Optimization Algorithms (2013, Wiley))
  while(parent < couples_nb*2){

    r = ((double) rand() / (RAND_MAX));

    for(int i = 0; i < N_moleculas; i++)
      if( soma_prob_new[i] < r && r < prob[i] + soma_prob_new[i] && is_parent[i] == 0){
        is_parent[i] = 1;
        parent_order[parent] = i;
        parent += 1;
      }
  }
}

void generate_children(vector<molecula*> pop, int *parent_order ){

    double r;

    srand((unsigned int)time(NULL));

    for(int i = 0; i < couples_nb; i++)

      for (int j = 0; j < children_per_couple; j++){

        double gene_prop = 1./(j+2);

        pop[parents_nb + i*children_per_couple + j] -> Mating(pop[ parent_order[2*i] ], pop[ parent_order[2*i + 1] ], gene_prop);

        r = ((double) rand() / (RAND_MAX));

        if ( r < mutation_prob)
          pop[parents_nb + i*children_per_couple + j] -> Mutate_1Atom();
        
      }
}


int main(){

  TCanvas *c1 = new TCanvas();
  auto gr = new TGraph();
  
  //population of molecules
  vector<molecula*> pop;

  for(int i = 0; i < N_moleculas; ++i)
    pop.push_back(new molecula(N_atomos, dim_caixa));
 
  // is_parent[i] == 0 : False: Molecule i is not a parent
  // is_parent[i] == 1 : True: Molecule i is not a parent
  bool *is_parent = new bool[N_moleculas];
  int *parent_order = new int[couples_nb*2];

  for(int iter = 0; iter < max_iter; iter++){

    for(int i = 0; i < N_moleculas; ++i)
      pop[i] -> Potencial();

    sort(pop.begin(), pop.end(), molecula::LessPot);

    gr -> AddPoint( iter, pop[0] -> Get_Pot());

    parent_probability( pop, is_parent, parent_order);
    generate_children( pop, parent_order);
  }

  delete[] is_parent;
  delete[] parent_order;

  pop.clear();

  c1 -> cd();
  gr -> Draw("AP");
  c1 -> SaveAs("evolution.pdf");

  return 0;
}
