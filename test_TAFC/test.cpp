#include <iostream>
#include <vector>
#include <fstream>

#include "TCanvas.h"
#include "TGraph.h"

#include "Molecule.h"

using namespace std;

int main(){
  unsigned short int N_mol = 30, N_atoms = 13, L_box = 10;
  int N_gens = 50000;
  float mute_rate = 0.1, surv_rate = 0.2;

  vector<Molecule*> pop;

  TCanvas* c1 = new TCanvas();
  TGraph* gr = new TGraph();

  for(int i=0; i<N_mol; ++i){
    pop.push_back(new Molecule(N_atoms, L_box, mute_rate));
    pop[i]->Fit();
  }

  sort(pop.begin(), pop.end(), Molecule::LessPot);

  for(int gen=1; gen<N_gens; ++gen){
    gr->AddPoint(gen, pop[0]->Get_Fit());
    
    for(int mol=surv_rate*N_mol; mol<N_mol; mol+=surv_rate*N_mol){
      for(int alive=0; alive<surv_rate*N_mol; ++alive)
	pop[mol+alive]->Set_Pos(pop[alive]->Get_Pos());
    }

    for(int mol=0; mol<N_mol; ++mol){
      pop[mol]->Mutate();
      pop[mol]->Fit();
    }

    sort(pop.begin(), pop.end(), Molecule::LessPot);

  }

  ofstream file("Best_Molecule.bs");

  double** best_pos = pop[0]->Get_Pos();
  for(int i=0; i<N_atoms; ++i){
    file << "atom C " << flush;
    for(int j=0; j<3; ++j){
      cout << best_pos[i][j] << " || " << flush;
      file << best_pos[i][j] << " " << flush;
    }
    cout << endl;
    file << endl;
  }

  //Check distances
  
  for(int i=0; i<N_atoms-1; ++i){
    for(int j=i+1; j<N_atoms; ++j){
      double distance = 0.;
      for(int k=0; k<3; ++k)
	distance += (best_pos[i][k] - best_pos[j][k])*(best_pos[i][k] - best_pos[j][k]);
      distance = sqrt(distance);
      cout << distance << endl;
    }
  }

  file << endl << "spec C 0.3 0.7" << endl
       << endl << "bonds C C 0.5 2.0 0.1 0.0"
       << endl << "bonds C H 0.4 1.0 0.1 1.0" << flush;
  
  gr->Draw("AP");
  c1->SaveAs("mutation.pdf");

  file.close();
  
  return 0;
}
