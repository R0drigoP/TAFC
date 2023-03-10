#include "molecula.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <ctime>

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
  for(int i = 0; i < N_moleculas; ++i)
    pop.push_back(new molecula(N_atomos, dim_caixa));

  //4 parents -> 2 couples
  //r = ((double) rand() / (RAND_MAX))
  // 3 children for each
  int parents_nb = int(survival_rate * N_moleculas + 0.5);
  //se for ímpar descartamos uma molecula como parent; há sempre bons partidos que acabam solteiros :(
  int couples_nb = int(parents_nb/2);
  cout << "Nr casais" << couples_nb << endl;
  int children_nb = N_moleculas - parents_nb;
  //temos de ver como vamos distribuir as children por casal
  int children_per_couple = children_nb / couples_nb;
  cout << children_per_couple << endl;


  TCanvas *c1 = new TCanvas();
  auto gr = new TGraph();

  double *prob = new double[N_moleculas];
  double soma_prob, soma_prob_aux;
  double *soma_prob_new = new double[N_moleculas];
  bool *is_parent = new bool[couples_nb*2];
  double r, s;

  srand((unsigned int)time(NULL));

  int nr_couples_formed = 0;

  int parent = 0;

  int *parent_order = new int[couples_nb*2];

  for(int iter = 0; iter < 1000; iter++){

    //cout<<"ITER "<<iter<<endl;

      for(int i = 0; i < N_moleculas; ++i)
        pop[i]->Potencial();

      sort(pop.begin(), pop.end(), molecula::LessPot);

      gr -> AddPoint( iter, pop[0] -> Get_Pot());
      //cout << "|| " << pop[0] -> Get_Pot() << endl;

      //imprimir
      if(iter == 0 || iter == (10-1) ){
      	cout << "ITER " << iter << endl;

	    for(int i = 0; i < N_moleculas; ++i)
	       cout << pop[i] -> Get_Pot() << endl;
      }

      //calculate probability for each molecule to be a parent
      double Emax = pop[N_moleculas-1] -> Get_Pot();
      double Emin = pop[0] -> Get_Pot();


      // is_parent[i] == 0 : False: Molecule i is not a parent
      // is_parent[i] == 1 : True: Molecule i is not a parent
      for(int i = 0; i < N_moleculas; i++) 
      	is_parent[i] = 0;

      soma_prob = 0;
      soma_prob_aux = 0;

      for(int i = 0; i < N_moleculas; i++){
      	soma_prob_new[i] = 0;
      }

      for(int i = 0; i < N_moleculas; i++){
      	//biblio: A Genetic Algorithm for Lennard-Jones Atomic Clusters C. BARR6N
      	//prob[i] = 1 + N_moleculas*pow(((Emax - pop[i] -> Get_Pot())/(Emax - Emin)),2);
      	prob[i] =  exp(-(pop[i]->Get_Pot()));

      	soma_prob += prob[i]; 
      }

      for(int i = 0; i < N_moleculas; i++){

      	prob[i] = prob[i]/soma_prob;

      	//cout << "Probabilidade molecula " << i << " : " << prob[i] << endl;

      	soma_prob_new[i] += soma_prob_aux; 

      	soma_prob_aux += prob[i];

      	//cout << "Soma prob:" << soma_prob_new[i] << endl;
      	
      }

  	//select couples of parents (de acordo c/ pagina 76/77 Evolutionary Optimization Algorithms)
    //!! Acho que os pais deviam ser replaced pelas suas children (Bibliografia: Evolutionary Optimization Algorithms (2013, Wiley))

  	while(parent < couples_nb*2){

  		r = ((double) rand() / (RAND_MAX));

  		for(int i = 0; i < N_moleculas; i++){

	  		if( soma_prob_new[i] < r && r < prob[i] + soma_prob_new[i] && is_parent[i] == 0){

	  			//cout << "[" << soma_prob_new[i] << ", "<< prob[i] + soma_prob_new[i] << "]" << endl;

	  			//cout << "random nr: " << r << "molecule nr: " << i << endl;
				//molecula i é parent 
				is_parent[i] = 1;

				parent_order[parent] = i;

				parent += 1;

			}
  		}
  	}


    for(int i = 0; i < couples_nb; i++){

        for (int j = 0; j < children_per_couple; j++){

          double gene_prop = 1./(j+2);

          pop[parents_nb + i*children_per_couple + j] -> Mating(pop[ parent_order[2*i] ], pop[ parent_order[2*i + 1] ], gene_prop);

          //cout << "Estou a morrer?" << i << j << endl;

          //cout << parents_nb+i*children_per_couple+j << endl;//cout<<couples_nb*i<<" "<<couples_nb*i+1<<endl;
        
        }
    }
  }

  delete[] prob;
  delete[] soma_prob_new;
  delete[] is_parent;
  delete[] parent_order;

  pop.clear();

  c1 -> cd();
  //gr->GetYaxis()->SetRangeUser(-1.,0.);
  gr -> Draw("AP");
  c1 -> SaveAs("evolution.pdf");

  cout << "ola" << endl;
  return 0;
}
