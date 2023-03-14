#include "molecula.h"

//Constructors and Destructor

molecula::molecula(int n_atomos, int dim_caixa, double mutation_rate) : N_atomos(n_atomos), Dim_caixa(dim_caixa), mute_rate(mutation_rate), f_value(0.) {
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

molecula::molecula(const molecula& mol) :
  molecula(mol.N_atomos, mol.Dim_caixa, mol.mute_rate) {;}

//Neste aqui é preciso passar nmr de atomos como argumento? Não se pode só usar um getter?
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

void molecula::Mating_Plano(molecula* mom, molecula* dad){

  cout << "Ola" << endl;
  //escolher parte dos atomos da mãe e parte dos átomos do pai
  double **pos_mom = mom -> Get_Pos();
  double **pos_dad = dad -> Get_Pos();

  double pos_CM[3];

  for(int i = 0; i < 3; i++)
    pos_CM[i] = 0;
  
  //posicao centro de massa 
  for(int i = 0; i < N_atomos; i++)
    for(int j = 0; j < 3; j++)
      pos_CM[j] += (pos_mom[i][j] + pos_dad[i][j])/2;

  for(int i = 0; i < 3; i++)
    pos_CM[i] = pos_CM[i]/N_atomos; 
  
  
  //cortar a caixa c/ plano c/ z cte.
  //Molecula filha resulta dos atomos da mãe que se encontram acima do plano + atomos do pai a baixo do plano
  int nr_atoms_mom = 0, nr_atoms_dad = 0, nr_atoms = 0;

  while( nr_atoms < N_atomos ){

    cout << "Ola" << endl;

    for (int i = 0; i < N_atomos; i++){

      if(pos_mom[i][2] > pos_CM[2]){
        for(int j = 0; j < 3; j++)
          posicoes[nr_atoms][j] = pos_mom[i][j];

        nr_atoms_mom += 1;
        nr_atoms += 1;
        cout << "Nr total atomos filho: "<<nr_atoms << endl;
        cout << "Nr atomos filho (mae): "<<nr_atoms_mom << endl;

      }

      if(nr_atoms == N_atomos )
        break;

      if(pos_dad[i][2] < pos_CM[2]){
        for(int j = 0; j < 3; j++)
          posicoes[nr_atoms][j] = pos_dad[i][j];

        nr_atoms_dad += 1;
        nr_atoms += 1;
        cout << "Nr total atomos filho: "<< nr_atoms << endl;
        cout << "Nr atomos filho (pai): "<<nr_atoms_dad << endl;
      }

      if(nr_atoms == N_atomos )
        break;

    if( i == (N_atomos-1) && nr_atoms < N_atomos){

      cout << "Ups!!!" << endl;

      if(nr_atoms == 0)
        cout << "Nao encontrei nada :( " << endl;

      for (int i = 0; i < N_atomos; i++){
        if(pos_mom[i][2] < pos_CM[2]){

          for(int j = 0; j < 3; j++)
            posicoes[nr_atoms][j] = pos_mom[i][j];

          nr_atoms_mom += 1;
          nr_atoms += 1;
          cout << "Nr total atomos filho: "<<nr_atoms << endl;
          cout << "Nr atomos filho (mae): "<<nr_atoms_mom << endl;

        }

        if(nr_atoms == N_atomos)
          break;

        if(pos_dad[i][2] > pos_CM[2]){

          for(int j = 0; j < 3; j++)
            posicoes[nr_atoms][j] = pos_dad[i][j];

          nr_atoms_dad += 1;
          nr_atoms += 1;
          cout << "Nr total atomos filho: "<< nr_atoms << endl;
          cout << "Nr atomos filho (pai): "<<nr_atoms_dad << endl;
        }

        if(nr_atoms == N_atomos )
          break;

    }
  }
  }
}
}


molecula::~molecula() {
    //apagar posicoes
    for(int i = 0; i < 3; ++i) 
      delete[] posicoes[i];
    
    delete[] posicoes;
}

//Operators
void molecula::operator=(const molecula& mol){
  for(int i=0; i<N_atomos; ++i){
    for(int j=0; j<3; ++j)
      posicoes[i][j] = mol.posicoes[i][j];
  }
}

//Setters

void molecula::Set_Pos(double** new_pos){
  for(int i=0; i<N_atomos; ++i){
    for(int j=0; j<3; ++j)
      posicoes[i][j] = new_pos[i][j];
  }
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

void molecula::OtherPotential(){
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
  
  f_value = f;
}


void molecula::Mutate(){
  gRandom = new TRandom3(0);
  double check_if_mute = gRandom->Uniform(0,1);
  
  if(check_if_mute < mute_rate){
    int atom_to_mutate = (int)gRandom->Uniform(0, N_atomos);
    for(int i=0; i<3; ++i){
      double mutation = gRandom->Uniform(-1,1)*0.01*Dim_caixa;
      posicoes[atom_to_mutate][i] += mutation;
      if(posicoes[atom_to_mutate][i] > Dim_caixa)
	posicoes[atom_to_mutate][i] -= 2*mutation;
    }
  }
}

void molecula::Mutate_1Atom(){
  gRandom = new TRandom3(0);
  int atom_to_mutate = int(gRandom->Uniform(0,N_atomos));
  
  for(int j=0; j<3; ++j){
    double mutation = gRandom->Uniform(-1,1)*0.01*Dim_caixa;
    posicoes[atom_to_mutate][j] += mutation;
    if(posicoes[atom_to_mutate][j] > Dim_caixa)
      posicoes[atom_to_mutate][j] -= 2*mutation;
    
  }
}
