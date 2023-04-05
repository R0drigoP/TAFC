#include "nr3.h"
#include "mins.h"
#include "mins_ndim.h"

#include "molecule.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include <ctime>
#include <fstream>
#include "global.h"
//

using namespace std;

//global variables
unsigned int N_molecules = 100, N_atoms = 38;
float L_box = 3., survival_rate = 0.7, mutation_prob = 0.01, sex_prob = 0.6;
float alpha = 0., m0 = 0.1;
unsigned int max_iter = 100000;
int N_times = 10;

unsigned int parents_nb = int(survival_rate * N_molecules);
//int couples_nb = int(parents_nb/2);
//int children_per_couple = ( N_moleculas - parents_nb) / couples_nb;

bool mating = 1;

unsigned int nb_of_calls = 0, nb_of_calls_mute = 0, nb_of_calls_mat = 0, nb_of_calls_mat_plano = 0;

double final_fit = 0.;

struct Funcd {
  
    Funcd(){}

    Doub operator() (VecDoub_I &x) {
      Doub f = 0.;
      for (int i = 0; i < N_atoms-1; ++i){
          for(int j = i+1; j < N_atoms; ++j){
              // Calculate radius
              Doub r = 0.;
              for(int k = 0; k < 3; ++k)
                  r += (x[i*3+k] - x[j*3+k])*(x[i*3+k] - x[j*3+k]);
              r = sqrt(r);

              f += 4*( pow(1./r, 12) - pow(1./r, 6) );
          }
      }
      return f;
    }

    void df(VecDoub_I &x, VecDoub_O &deriv) {

        for (int i = 0; i < N_atoms-1; ++i){
            for(int j = i+1; j < N_atoms; ++j){
                // Calculate radius
                Doub r = 0.;
                for(int k = 0; k < 3; ++k)
                    r += (x[i*3+k] - x[j*3+k])*(x[i*3+k] - x[j*3+k]);
                r = sqrt(r);

                // Calculate force and add to deriv
                Doub force = -24*( 2*pow(1./r, 13) - pow(1./r, 7) );
                for(int k = 0; k < 3; ++k){
                    deriv[i*3+k] += force*(x[i*3+k] - x[j*3+k]);
                    deriv[j*3+k] += force*(x[j*3+k] - x[i*3+k]);
                }
            }
        }
    }
};

//probability of each molecule to be a parent
int main(){
  //cout << "ola1 " << endl;

  if(mating==1 &&  N_molecules * survival_rate < 2.){
    cout<<"To have sexual reprodution at least 2 molecules must survive each gen..."<<endl;
    cout<<"Increase your population or the survival probability"<<endl;
    return 1;
  }

  if(survival_rate < 0. || survival_rate > 1. || parents_nb >= N_molecules){
    cout << "Check your survival_rate" << endl;
    return 1;
  }

  double min_pot = 0.;
  int min_nb_calls = 0;
  for(int n=0; n<N_times; ++n){
    nb_of_calls = 0;
    double **positions;
    positions = new double*[N_atoms];
    for (int i = 0; i < N_atoms; i++) 
      positions[i] = new double[3];

    int* flag = new int[N_molecules];
    for (int i = 0; i < N_molecules; i++) 
      flag[i] = 0;


    TCanvas *c1 = new TCanvas();
    auto gr = new TGraph();

    ofstream text_file("best_molecule.bs");
    ofstream movie_file("best_molecule.mv");
  
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
    
    for(int iter = 0; iter < max_iter; iter++){

      //cout << "ITER NR " << iter << endl;
      /*
	for(int i = 0; i < N_moleculas; i++){
	cout<<"Molecula"<<i<<endl;
	positions = pop[i] -> Get_Pos();
        cout << "Atomo " << 0 << " : " <<positions[0][0] << ", " << positions[0][1] << ", " << positions[0][2] << endl;
        cout << "Atomo " << 1 << " : " <<positions[1][0] << ", " << positions[1][1] << ", " << positions[1][2] << endl;
        cout << "Atomo " << 2 << " : " <<positions[2][0] << ", " << positions[2][1] << ", " << positions[2][2] << endl;

	}*/

      //cout<<"got potential"<<endl;

      sort(pop.begin(), pop.end(), molecule::LessPot);

      if(iter == 0){
	double** best_pos = pop[0]->Get_Pos();
	for(int i = 0; i < N_atoms; ++i){
	  text_file << "atom C " << flush;
	  for(int j = 0; j < 3; ++j)
	    text_file << best_pos[i][j] << " " << flush;
	  text_file << endl;
	}
      }
      else if(iter % 1000 == 0){
	movie_file << "frame" << endl;
	double** best_pos = pop[0]->Get_Pos();
	for(int i = 0; i < N_atoms; ++i){
	  for(int j = 0; j < 3; ++j)
	    movie_file << best_pos[i][j] << " ";
	}
	movie_file << endl << endl;

	//cout << "ITER NR " << iter << endl;

	//cout << "Pot so far " << pop[0]->Get_Fit() << endl;

	//cout<<"nb of mutations "<<nb_of_calls_mute<<endl;
	//cout<<"nb of mat plano "<<nb_of_calls_mat_plano<<endl;
	//cout<<"nb of mat       "<<nb_of_calls_mat<<endl;
	//cout<<"pop x iter      "<<N_molecules*(iter+1)<<endl;
	//cout<<"total reproduti "<<nb_of_calls_mute+nb_of_calls_mat_plano+nb_of_calls_mat<<endl;
	//cout<<"nb of func calss"<<nb_of_calls<<endl;

      }
    
      //gr -> AddPoint( iter, pop[0] -> Get_Fit());
 
      //matar os mais fracos -> fazer copias da melhor pop (se calhar atribuir alguma aleatoriadade a este processo)
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
    
	if( mating == 1){
	  //setting flag to 0 for parents
	  for(int i = 0 ; i < parents_nb; i++)
	    flag[i] = 0;
      
	  //sexual reproduction
	  TRandom3* gRandom = new TRandom3(0); 
#pragma omp for
	  for(int i = parents_nb ; i <N_molecules; i++)
	    flag[i] = pop[i]->generate_children3(pop, gRandom);
      
	  delete gRandom;
	}
    
	//assexual reproduction
	TRandom3* gRandom = new TRandom3(0); 
    
#pragma omp for
	for(int mol = 0; mol < N_molecules; mol++){
	  pop[mol] -> Mutate(iter, m0, alpha, flag[mol], gRandom);
	}
      
	delete gRandom;
      }
    
    
      //cout << "Pot: " << pop[0] -> Get_Fit() << endl;

      if(iter == max_iter-1){
	positions = pop[0] -> Get_Pos();
	final_fit =  pop[0]-> Get_Fit();
	//cout << "Final Pot " << final_fit << endl;
      
	//for(int i = 0; i < N_atomos; i++)
	//cout << "Atomo " << i << " : " << positions[i][0] << ", " << positions[i][1] << ", " << positions[i][2] << endl;
      }
    
    }

    positions = (pop[0] -> Get_Pos());
    
    //so podemos fazer isto quando encontrarmos o max global
     for(int i = 0; i < N_atoms; i++ )
      for(int j = 0; j < 3; j++){
        //cout << "ola" << i+5 << endl;
        p[i*3+j] = positions[i][j];
      }
    
  
    pop[0] -> Fit();

    //cout << "Pot: " << pop[0] -> Get_Fit() << endl;

    p = frprmn.minimize(p);

    for(int i = 0; i < N_atoms; i++ )
      for(int j = 0; j < 3; j++)
	positions[i][j] = p[i*3+j];
    
    pop[0] -> Set_Pos(positions);

    pop[0] -> Fit();

    if(pop[0]->Get_Fit() < min_pot){
      min_pot = pop[0]->Get_Fit();
      min_nb_calls = nb_of_calls;
    }
	
    //cout << "Pot: " << pop[0] -> Get_Fit() << endl;
    //delete[] is_parent;
    //delete[] parent_order;
  
    for(int i = 0; i < 3; ++i) 
      delete[] positions[i];
  
    delete[] positions;
    delete[] flag;
    //Forma correta de destruir o vetor mas dá seg fault
    //for(vector<molecule*>::iterator it = pop.begin(); it != pop.end(); ++it)
    //delete *it;
    
    pop.clear();
  
  
    double atom_size = 0.1/L_box;
    
    text_file << endl << "spec C 0.1 Red" << endl
	      << endl << "bonds C C 0.5 1.5 0.01 0.0"
	      << endl << "bonds C H 0.4 1.0 0.01 1.0"
	      << endl << "scale 100"
	      << endl << "inc 5.0" << flush;
    
    //c1 -> cd();
    //gr->GetHistogram()->SetMaximum(-0.1*final_fit);
    //gr->GetHistogram()->SetMinimum(1.01*final_fit);
    
    //gr -> Draw("AP");
    //c1 -> SaveAs("evolution.pdf");
    
    text_file.close();
    movie_file.close();
    
  }
  
  cout << "Min value is " << min_pot << endl;
  
  cout<<"Total pot calls: "<<min_nb_calls<<endl;
  
  return 0;
}
