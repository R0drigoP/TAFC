#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "molecula.h"
#include "TCanvas.h"
#include "TGraph.h"

using namespace std;

int main()
{
  int N_moleculas = 30;
  int N_atomos = 13;
  int dim_caixa = 10;
  int N_gen = 10;
  double mutation_rate = 0.01; //5%
  double survival_rate = 0.3; //40%

  vector<double> Pot_t(N_gen); //Potencial over generations
  
  //Inicializar população como vector de ponteiros de objetos
  vector<molecula*> pop;

  //Initial population and calculate their potential
  for(int i=0; i<N_moleculas; ++i){
    pop.push_back(new molecula(N_atomos, dim_caixa, mutation_rate));
    pop[i]->Potencial();
  }
  

  for(int gen=0; gen<N_gen; ++gen){
    sort(pop.begin(), pop.end(), molecula::LessPot);

    
    cout << "\nGeneration number" << gen + 1 << endl;
    //cout << pop[0]->Get_Pot() << endl;
    
    for(int mol=0; mol<N_moleculas; ++mol)
      cout << pop[mol]->Get_Pot() << " ||  " << flush;

    Pot_t[gen] = pop[0]->Get_Pot();

    //Nova gen formada
    for(int mol=survival_rate*N_moleculas; mol<N_moleculas; mol+=survival_rate*N_moleculas){
      for(int alive=0; alive<survival_rate*N_moleculas; ++alive)
	pop[mol+alive] = pop[alive];
    }

    //Introduzir mutação
    for(int mol=0; mol<N_moleculas; ++mol){
      pop[mol]->Mutate();
      pop[mol]->Potencial();
    }
  }

  TCanvas *c1 = new TCanvas();
  auto gr = new TGraph();
  
  for(int gen=0; gen<N_gen; ++gen){
    gr->AddPoint(gen+1,Pot_t[gen]);
    //cout << Pot_t[gen] << endl;
  }

  gr->Draw("AP");
  c1->SaveAs("mutation.pdf");
  //Fazer o clear das cenas
    

  return 0;
}