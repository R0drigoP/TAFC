#include "molecula.h"

//Constructors and Destructor

molecula::molecula(int n_atomos, int dim_caixa) : N_atomos(n_atomos), Dim_caixa(dim_caixa){
  posicoes = new double*[n_atomos];
  
  for (int i = 0; i < n_atomos; i++) {
    posicoes[i] = new double[3];
  }
  
  gRandom = new TRandom3(0);
  for (int i = 0; i < N_atomos; i++){
    for (int j = 0; j < 3; j++){
    posicoes[i][j] = gRandom -> Uniform(-dim_caixa, dim_caixa); //gera posicao entre -Dim e Dim
    //cout << "Pos: x:"<< posicoes[j][0] << " - y:" << posicoes[j][1] << " - z:"  << posicoes[j][2] << endl;
    }
  }
}

molecula::molecula(molecula* mom, molecula* dad, double gene_prop, int n_atomos): N_atomos(n_atomos) {
  double **pos_mom = mom->Get_Pos();
  double **pos_dad = dad->Get_Pos();

  for (int i = 0; i < N_atomos; i++){
    for (int j = 0; j < 3; j++){
      posicoes[i][j] = gene_prop * pos_mom[i][j] + (1 - gene_prop)* pos_dad[i][j];
      }
    }
}

void molecula::Mating(molecula* mom, molecula* dad, double gene_prop){
  double **pos_mom = mom->Get_Pos();
  double **pos_dad = dad->Get_Pos();

  for (int i = 0; i < N_atomos; i++){
    for (int j = 0; j < 3; j++){
      posicoes[i][j] = gene_prop * pos_mom[i][j] + (1 - gene_prop)* pos_dad[i][j];
      }
    }
}



molecula::~molecula() {
    //apagar posicoes
    for(int i = 0; i < 3; ++i) { delete[] posicoes[i]; }
    
    delete[] posicoes;
}


//calcula o lennard_jones para cada molecula
double molecula::Potencial(){
  float sigma = 1; //(Angstroms)
  float epsilon = 1; //(kJ/mol)
  
  float d_aux = 0;

  f_value = 0.;
  
  //nao sei se esta certo, confirmar com os outros
  for (int j = 0; j < N_atomos; j++){
    
    for (int k = 0; k < N_atomos; k++){

      if( j < k ){
        d_aux = sqrt( pow((posicoes[j][0] - posicoes[k][0]), 2) + pow((posicoes[j][1] -posicoes[k][1]), 2) + pow((posicoes[j][2] - posicoes[k][2]), 2));

        f_value += 4*epsilon*(pow((sigma/d_aux), 12) - pow((sigma/d_aux), 6));		
      }
    }
  }
  
  //cout << "f_value: "<< f_value << endl;
  
  return f_value;
}

double molecula::OtherPotential(){
  float sigma = 1; //(Angstroms)
  float epsilon = 1; //(kJ/mol)
  double f = 0.;

  for(int i=0; i<N_atomos-1; ++i){
    for(int j=i+1; j<N_atomos; ++j){
      //Calculate radius
      double r = 0.;
      for(int k=0; k<3; ++k)
	r += (posicoes[i][k] - posicoes[j][k])*(posicoes[i][k] - posicoes[j][k]);
      r = sqrt(r);

      //Calculate Potential and sum to f_value
      f += 4*epsilon*(pow((sigma/r), 12) - pow((sigma/r), 6));
    }
  }
  
  return f;
}


















