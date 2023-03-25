#include "mating_func.h"
#include "global.h"
#include <ctime>

void parent_probability(vector<molecula*> pop, bool *is_parent, int *parent_order ){

  //cout<<"calculating parent_probability"<<endl;
  
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

  //cout<<"pos os vecs a 0"<<endl;

  for(int i = 0; i < N_moleculas; i++){
    prob[i] =  exp( -(pop[i] -> Get_Pot()) );
    //cout << pop[i] -> Get_Pot()<<endl;
    //cout<<exp( -(pop[i] -> Get_Pot()) )<<endl;
    soma_prob += prob[i]; 
  }

  //cout<<"calculou probs"<<endl;

  for(int i = 0; i < N_moleculas; i++){
    prob[i] = prob[i]/soma_prob;
    soma_prob_new[i] += soma_prob_aux; 
    soma_prob_aux += prob[i];
  }

  //for(int i = 0; i < N_moleculas; i++){
  //  cout<<prob[i]<<endl;
  //}



  while(parent < couples_nb*2){

    //cout<<"a"<<endl;

    r = ((double) rand() / (RAND_MAX));

    for(int i = 0; i < N_moleculas; i++){

      //cout<< soma_prob_new[i] << "  "<< r<<"  " << soma_prob_new[i]<<endl;

      if( soma_prob_new[i] < r && r < prob[i] + soma_prob_new[i] && is_parent[i] == 0){

        is_parent[i] = 1;
        parent_order[parent] = i;
        parent += 1;
        //cout<<parent<<endl;
      }
    }
  }
  //cout<<"calculated parent_probability"<<endl;

  //delete[] prob;
  //delete[] soma_prob_new;
}

void generate_children(vector<molecula*> pop, int *parent_order ){

    double r;

    srand((unsigned int)time(NULL));

    for(int i = 0; i < couples_nb; i++)

      for (int j = 0; j < children_per_couple; j++){
        
        double gene_prop = 1./(j+2);

        pop[parents_nb + i*children_per_couple + j] -> Mating_Plano3(pop[ parent_order[2*i] ], pop[ parent_order[2*i + 1] ]);
        //pop[parents_nb + i*children_per_couple + j] -> Mating_Plano3(pop[ parent_order[2*i] ], pop[ parent_order[2*i + 1] ]);


        //pop[parents_nb + i*children_per_couple + j] -> Mutate();
        
      }
}

void generate_children2(vector<molecula*> pop){

    double r;

    srand((unsigned int)time(NULL));

    for(int i = 0; i < couples_nb; i++)

      for (int j = 0; j < children_per_couple; j++){
        
        double gene_prop = 1./(j+2);

        pop[parents_nb + i*children_per_couple + j] -> Mating_Plano3(pop[ couples_nb*i ], pop[ couples_nb*i+1 ]);
        //pop[parents_nb + i*children_per_couple + j] -> Mating_Plano3(pop[ parent_order[2*i] ], pop[ parent_order[2*i + 1] ]);


        //pop[parents_nb + i*children_per_couple + j] -> Mutate();
        
      }
}

