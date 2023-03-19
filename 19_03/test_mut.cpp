#include <cstdio>
#include <iostream>
#include <algorithm>
#include <vector>

#include "molecula.h"
#include "TCanvas.h"
#include "TGraph.h"

using namespace std;

int main(){
  int N_mol = 10, N_atoms = 10, L_box = 10, N_gens = 500;
  float mute_rate = 0.001, surv_rate = 0.5;
  vector<double> best_pot(N_gens); //vector q guarda melhores potenciais
  
  
  vector<molecula*> pop(N_mol);

  TCanvas* c1 = new TCanvas();
  auto gr = new TGraph();
  
  for(int i=0; i<N_mol; ++i){
    pop[i] = new molecula(N_atoms, L_box, mute_rate);
    //pop[i].Potencial();
  }
  
  //cout << pop[0]->Potencial() << endl;
  
  
  for(int gen=0; gen<N_gens; ++gen){
    sort(pop.begin(), pop.end(), molecula::LessPot);

    for(int i=0; i<N_atoms; ++i)
      cout << pop[i]->Get_Pot() << " || " << flush;
    cout << endl;
    
    best_pot[gen] = pop[0]->Get_Pot();
    gr->AddPoint(gen+1,best_pot[gen]);
    
    for(int mol=surv_rate*N_mol; mol<N_mol; mol+=surv_rate*N_mol){
      //cout << mol << endl;
      for(int alive=0; alive<surv_rate*N_mol; ++alive){
	double x = 1.;
	pop[mol+alive] = pop[alive];
      }
    }
    
    for(int mol=0; mol<N_mol; ++mol){
      pop[mol]->Mutate_1Atom();
      pop[mol]->Potencial();
    }
    
    //cout << best_pot[gen] << endl;
  }

  gr->Draw("AP");
  c1->SaveAs("test_mut.pdf");
  //cout << best_pot[N_gens-1] << endl;
  //cout << "End of cycle" << endl;
  
  return 0;
}
