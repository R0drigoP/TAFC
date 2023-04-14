

#include "molecule.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include <ctime>
#include <fstream>
#include "global.h"

//#include <omp.h>


//

using namespace std;

//global variables


unsigned int N_molecules = 1, N_atoms = 38;
float L_box = 2., survival_rate = 0.95, mutation_prob = 0.01, sex_prob = 0.5;
float alpha = 0., m0 = 0.5;
unsigned int max_iter = 1000000;


unsigned int parents_nb = int(survival_rate * N_molecules);
//int couples_nb = int(parents_nb/2);
//int children_per_couple = ( N_moleculas - parents_nb) / couples_nb;

bool mating = 0;

unsigned int nb_of_calls = 0, nb_of_calls_mute = 0, nb_of_calls_mat = 0, nb_of_calls_mat_plano = 0;

double final_fit = 0.;

/*
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
};*/

//probability of each molecule to be a parent
int main(){

  //double t0 = omp_get_wtime();
  /*
  if(mating==1 &&  N_molecules * survival_rate < 2.){
    cout<<"To have sexual reprodution at least 2 molecules must survive each gen..."<<endl;
    cout<<"Increase your population or the survival probability"<<endl;
    return 1;
  }

  if(survival_rate < 0. || survival_rate > 1. || parents_nb >= N_molecules){
    cout << "Check your survival_rate" << endl;
    return 1;
  }*/

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

  
    double** best_pos = pop[0]->Get_Pos();
    for(int i = 0; i < N_atoms; ++i){
      text_file << "atom C " << flush;
      for(int j = 0; j < 3; ++j)
        text_file << best_pos[i][j] << " " << flush;
      text_file << endl;
    }

    for (int i=1; i<N_molecules;i++){
      movie_file << "frame" << endl;
      best_pos = pop[i]->Get_Pos();
      for(int i = 0; i < N_atoms; ++i){
        for(int j = 0; j < 3; ++j)
          movie_file << best_pos[i][j] << " ";
      }
      movie_file << endl << endl;

    }
  

  
  text_file << endl << "spec C 0.1 Red" << endl
	    << endl << "bonds C C 0.3 1.0 0.01 0.0"
	    << endl << "bonds C H 0.4 1.0 0.01 1.0"
	    << endl << "scale 100"
	    << endl << "inc 5.0" <<endl<< flush;
  /*

  c1 -> cd();
  gr->GetHistogram()->SetMaximum(-0.1*final_fit);
  gr->GetHistogram()->SetMinimum(1.01*final_fit);

  gr -> Draw("AP");
  c1 -> SaveAs("evolution.pdf");*/
  
  text_file.close();
  movie_file.close();


  return 0;
}
