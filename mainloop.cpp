#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "molecula.h"

#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

using namespace std;

int main()
{
  int N_moleculas = 10;
  int N_atomos = 10;
  int dim_caixa = 10;
  double survival_rate = 0.4;
  
  //Inicializar população como vector de ponteiros de objetos
  vector<molecula*> pop;

  //Preencher o vector
  for(int i=0; i<N_moleculas; ++i)
    pop.push_back(new molecula(N_atomos, dim_caixa));

  

  //4 parents -> 2 couples
  // 3 children for each
  int parents_nb = int(survival_rate * N_moleculas + 0.5);
  int couples_nb = parents_nb/2;
  int children_nb = N_moleculas - parents_nb;
  int children_per_couple = children_nb / couples_nb;


 
    //_________
  TCanvas *c1 = new TCanvas();
  auto gr = new TGraph();


  for(int iter=0; iter<2000; iter++){

    //cout<<"ITER "<<iter<<endl;

      for(int i=0; i<N_moleculas; ++i)
        pop[i]->Potencial();

      sort(pop.begin(), pop.end(), molecula::LessPot);

      gr->AddPoint(iter,pop[0]->Get_Pot());
      cout<<"|| "<<pop[0]->Get_Pot()<<endl;

      //imprimir
      if(iter==0 || iter==(2000-1)){
        cout<<"ITER "<<iter<<endl;
        for(int i=0; i<N_moleculas; ++i)
          cout << pop[i]->Get_Pot() << endl;
      }

      //parent survival probability
      //double Emax = pop[N_moleculas-1];
      //double Emin = pop[0];



      for(int i = 0; i<couples_nb;i++){

        for (int j = 0; j < children_per_couple; j++){

          double gene_prop = 1./(j+2);
          pop[parents_nb+i*children_per_couple+j]->Mating(pop[couples_nb*i],pop[couples_nb*i+1], gene_prop);

          //cout << parents_nb+i*children_per_couple+j << endl;//cout<<couples_nb*i<<" "<<couples_nb*i+1<<endl;
        }

      }  

  }

  c1->cd();
  //gr->GetYaxis()->SetRangeUser(-1.,0.);
  gr->Draw("AP");
  c1->SaveAs("evolution.pdf");
  return 0;
}
