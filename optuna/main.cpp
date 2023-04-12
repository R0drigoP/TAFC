
#include "molecule.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include <ctime>
#include <fstream>


//#include <omp.h>


//

using namespace std;

//global variables
unsigned int N_atoms = 13;




unsigned int parents_nb = int(survival_rate * N_molecules);
//int couples_nb = int(parents_nb/2);
//int children_per_couple = ( N_moleculas - parents_nb) / couples_nb;

unsigned int nb_of_calls = 0, nb_of_calls_mute = 0, nb_of_calls_mat = 0, nb_of_calls_mat_plano = 0;

double final_fit = 0.;


//probability of each molecule to be a parent
int main(){

  //double t0 = omp_get_wtime();

  if(mating==1 &&  N_molecules * survival_rate < 2.){
    cout<<"To have sexual reprodution at least 2 molecules must survive each gen..."<<endl;
    cout<<"Increase your population or the survival probability"<<endl;
    return 1;
  }

  if(survival_rate < 0. || survival_rate > 1. || parents_nb >= N_molecules){
    cout << "Check your survival_rate" << endl;
    return 1;
  }

  double **positions;
  positions = new double*[N_atoms];
  for (int i = 0; i < N_atoms; i++) 
    positions[i] = new double[3];

  int* flag = new int[N_molecules];
  for (int i = 0; i < N_molecules; i++) 
    flag[i] = 0;

  ofstream opt_file("opt.txt");
  ofstream opt_calls("calls.txt");

  
  //population of molecules
  vector<molecule*> pop(N_molecules);

  for(int i = 0; i < N_molecules; ++i)
    pop[i] = new molecule(N_atoms, L_box, mutation_prob);
 
  // is_parent[i] == 0 : False: Molecule i is not a parent
  // is_parent[i] == 1 : True: Molecule i is not a parent
  //bool *is_parent = new bool[N_moleculas];
  //int *parent_order = new int[couples_nb*2];

  Funcd funcd;
  Frprmn<Funcd> frprmn(funcd);
  VecDoub p(3*N_atoms);

  //variables for early stopping
  double best = 0.; //so far best potential
  int iter_stop = 0; //iterations with the same best so far
  int acceptance = 500; //nb of iterations allowed in the same best

  for(int iter = 0; iter < max_iter; iter++){


    //REPRODUCTION 

    #pragma omp parallel
    {
      //---sexual reproduction
      if( mating == 1){

       
        for(int i = 0 ; i < parents_nb; i++)
          flag[i] = 0; //setting flag to 0 for parents

        TRandom3* gRandom = new TRandom3(0); 


        #pragma omp for
        for(int i = parents_nb ; i <N_molecules; i++)
          flag[i] = pop[i]->generate_children3(pop, gRandom);

        delete gRandom;
      }

      //---assexual reproduction

      TRandom3* gRandom = new TRandom3(0); 
      #pragma omp for
      for(int mol = 0; mol < N_molecules; mol++)
        pop[mol] -> Mutate(iter, m0, alpha, flag[mol], gRandom);
      
      delete gRandom;

    }

    //sort population
    sort(pop.begin(), pop.end(), molecule::LessPot);

    //matar os mais fracos e fazer copias da melhor pop (se calhar atribuir alguma aleatoriadade a este processo)
    #pragma omp parallel
    {
      #pragma omp for
      for(int mol = static_cast<int>(survival_rate*N_molecules); mol < N_molecules; mol += static_cast<int>(survival_rate*N_molecules)){
        int alive = 0;
        while(alive < survival_rate*N_molecules && (mol+alive)<N_molecules){
        //cout<<mol<<" "<<alive<<endl;                          
          pop[mol+alive] -> Set_Pos(pop[alive] -> Get_Pos());
          pop[mol+alive] -> Set_Fit(pop[alive] -> Get_Fit());

          ++alive;
        }
      }
    }

    

    //prints
    /*
    if(iter%100 == 0)
      cout << "ITER NR " << iter << " Pot: " << pop[0] -> Get_Fit() << endl;
    if(iter == max_iter-1){
      positions = pop[0] -> Get_Pos();
      final_fit =  pop[0]-> Get_Fit();
      cout << "Final Pot " << final_fit << endl;
    }*/

    //early stopping
    double current = pop[0]->Get_Fit();
    if (iter==0)
      best = current;
    if ( current < best){
      iter_stop=0;
      best = current;
    }
    else{
      iter_stop ++;
      if(iter_stop == acceptance  )
        break;
    }
    
  }//closing iterations loop

  positions = (pop[0] -> Get_Pos());
  for(int i = 0; i < N_atoms; i++ )
    for(int j = 0; j < 3; j++)
      p[i*3+j] = positions[i][j];

  
  pop[0] -> Fit();

  cout << "Pot: " << pop[0] -> Get_Fit() << endl;

  /*p = frprmn.minimize(p);

  for(int i = 0; i < N_atoms; i++ )
    for(int j = 0; j < 3; j++)
      positions[i][j] = p[i*3+j];

  pop[0] -> Set_Pos(positions);

  pop[0] -> Fit();

  cout << "Pot: " << pop[0] -> Get_Fit() << endl;*/

  opt_file<<pop[0] -> Get_Fit()<<flush;
  opt_file.close();

  opt_calls<<nb_of_calls<<flush;
  opt_calls.close();
  
  for(int i = 0; i < 3; ++i) 
    delete[] positions[i];
  
  delete[] positions;
  delete[] flag;
  //Forma correta de destruir o vetor mas dá seg fault
  //for(vector<molecule*>::iterator it = pop.begin(); it != pop.end(); ++it)
    //delete *it;
  
  pop.clear();
  
  double atom_size = 0.1/L_box;
  

  
  


  cout<<"Total pot calls: "<<nb_of_calls<<endl;

  /*
  double t1 = omp_get_wtime();
  printf("\n");
  printf("Number of threads   =  %i\n", omp_get_max_threads());
  printf("Computation time    =  %f ms\n", (t1-t0) * 1000);  */
  return 0;
}
