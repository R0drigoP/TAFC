#include "mating_func.h"
#include "global.h"
#include <ctime>

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

        pop[parents_nb + i*children_per_couple + j] -> Mating_Plano(pop[ parent_order[2*i] ], pop[ parent_order[2*i + 1] ]);

        r = ((double) rand() / (RAND_MAX));

        if ( r < mutation_prob)
          pop[parents_nb + i*children_per_couple + j] -> Mutate_1Atom();
        
      }
}