#include "molecule.h"


//Constructors and Destructor

molecule::molecule(unsigned int n_atoms, float l_box, float mute_prob) : N_atoms(n_atoms), L_box(l_box), mutation_prob(mute_prob) {
  positions = new double*[n_atoms];
  
  for (int i = 0; i < n_atoms; i++)
    positions[i] = new double[3];
  
  gRandom = new TRandom3(0);
  
  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++)
    positions[i][j] = gRandom -> Uniform(0., l_box);
  
}
  delete gRandom;

  this->Fit();
}

//Neste aqui é preciso passar nmr de atomos como argumento? Não se pode só usar um getter?
molecule::molecule(molecule* mom, molecule* dad, double gene_prop, unsigned int n_atoms): N_atoms(n_atoms) {
  double **pos_mom = mom->Get_Pos();
  double **pos_dad = dad->Get_Pos();
  
  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++){
      positions[i][j] = gene_prop * pos_mom[i][j] + (1 - gene_prop)* pos_dad[i][j];
    }
  }
}

molecule::~molecule() {
  for(int i = 0; i < 3; ++i){ 
    delete[] positions[i];
  }
  
  delete[] positions;
}

//Private method that calculates fitness
void molecule::Fit(){
  double f = 0.;

  ++nb_of_calls;

  for(int i = 0; i < N_atoms-1; ++i){
    for(int j = i+1; j < N_atoms; ++j){
      //Calculate radius
      double r = 0.;
      for(int k = 0; k < 3; ++k)
        r += (positions[i][k] - positions[j][k])*(positions[i][k] - positions[j][k]);
      r = sqrt(r);

      //Calculate Potential and sum to f_value
      f += 4*( pow(1./r, 12) - pow(1./r, 6) );
    }
  }
  
  fitness = f;
}


//Setters

void molecule::Set_Pos(double** new_pos){
  for(int i=0; i<N_atoms; ++i){
    for(int j=0; j<3; ++j)
      positions[i][j] = new_pos[i][j];
  }
}

void molecule::Set_Fit(double fit){
  fitness = fit;
}

//Mutation and Reproduction

void molecule::Mutate(unsigned int iter, float m0, float alpha, int flag, TRandom3* gRandom){
  //gRandom = new TRandom3(0);
  double check_if_mute = gRandom->Uniform(0,1);
  
  if(check_if_mute < mutation_prob){
    nb_of_calls_mute ++;
    int atom_to_mutate = (int)gRandom->Uniform(0, N_atoms);
    for(int i=0; i<3; ++i){
      
      double x0 = gRandom->Uniform(-1,1)*m0*L_box;

      double mutation = x0/(1+alpha*iter);

      positions[atom_to_mutate][i] += mutation;
      if(positions[atom_to_mutate][i] < 0. || positions[atom_to_mutate][i] > L_box)
        positions[atom_to_mutate][i] -= 2*mutation;
    }
    this->Fit();
  }
  //se nao tiver feito mutacao calcula na mesma o potencial para os que sofreram reproducao sexuada
  else if(flag==1)
    this->Fit();

  //delete gRandom;
}


int molecule::generate_children3(vector<molecule*> pop, TRandom3* gRandom){

  //gRandom = new TRandom3(0);

  double check_if_sex = gRandom->Uniform(0,1);

  if(check_if_sex < sex_prob ){

    //choose random parents from the surviving population
    int mom_index = (int)gRandom->Uniform(0, parents_nb);
    int dad_index = (int)gRandom->Uniform(0, parents_nb);

    //and make sure they are different
    while(dad_index==mom_index)
      dad_index = (int)gRandom->Uniform(0, parents_nb);

    double mating_type = gRandom->Uniform(0,1);

    //call reproduction method
    //maybe add some new ones later...
    if(mating_type<0.9){
      nb_of_calls_mat_plano ++;
      this->Mating_Plano3(pop[mom_index],pop[dad_index], gRandom);
    }

    else{
      nb_of_calls_mat ++;
      this->Mating(pop[mom_index],pop[dad_index],gRandom->Uniform(0,1));
    }
    
    //delete gRandom;
    return 1;
  }
  
  //delete gRandom;   
  return 0;
  
}


void molecule::Mating(molecule* mom, molecule* dad, double gene_prop){
  double **pos_mom = mom->Get_Pos();
  double **pos_dad = dad->Get_Pos();
  
  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++){
      positions[i][j] = gene_prop * pos_mom[i][j] + (1 - gene_prop)* pos_dad[i][j];
    }
  }
}


void molecule::Mating_Plano3(molecule* mom, molecule* dad, TRandom3* gRandom){

  double **pos_mom = mom -> Get_Pos();
  double **pos_dad = dad -> Get_Pos();
  double pos_CM = 0;

  int flag = 0;
  int nr_atoms = 0;

  //gRandom = new TRandom3(0);
  int dir = (int)gRandom -> Uniform(0,3);

  
  //posicao centro de massa 
  for(int i = 0; i < N_atoms; i++)
    pos_CM += (pos_mom[i][dir] + pos_dad[i][dir])/2;
  pos_CM = pos_CM/N_atoms; 

  /*

  cout << "Mom Atom 1 " <<  " : " <<pos_mom[0][0] << ", " << pos_mom[0][1] << ", " << pos_mom[0][2] << endl;
  cout << "Mom Atom 2 " << " : " <<pos_mom[1][0] << ", " << pos_mom[1][1] << ", " << pos_mom[1][2] << endl;
  cout << "Mom Atom 3 " << " : " <<pos_mom[2][0] << ", " << pos_mom[2][1] << ", " << pos_mom[2][2] << endl;

  cout << "dad Atom 1 " <<  " : " <<pos_dad[0][0] << ", " << pos_dad[0][1] << ", " << pos_dad[0][2] << endl;
  cout << "dad Atom 2 " << " : " <<pos_dad[1][0] << ", " << pos_dad[1][1] << ", " << pos_dad[1][2] << endl;
  cout << "dad Atom 3 " <<  " : " <<pos_dad[2][0] << ", " << pos_dad[2][1] << ", " << pos_dad[2][2] << endl;

  cout<<"z_cm "<<pos_CM[2]<<endl;*/


  //choose every mom atom above CM
  for (int i = 0; i < N_atoms; i++){


    if(pos_mom[i][dir] > pos_CM){
      for (int k = 0; k < nr_atoms; k++){
        flag = 0;
        for (int j = 0; j < 3; j++){
          if (positions[k][j]==pos_mom[i][j])
            flag++;
        }
        if(flag == 3)
          break; 
        }
        if (flag != 3){
          for(int j = 0; j < 3; j++)
            positions[nr_atoms][j] = pos_mom[i][j];
          nr_atoms ++;
          //cout<<"mom  above"<<endl;
        }
    } 
  }

  //choose every dad atom bellow CM
  for (int i = 0; i < N_atoms; i++){
    if(nr_atoms== N_atoms)
      break;
    else if(pos_dad[i][dir] < pos_CM){
      for(int j = 0; j < 3; j++)
        positions[nr_atoms][j] = pos_dad[i][j];
      nr_atoms ++;
      //cout<<"dad bellow"<<endl;
    }
  }

  if (nr_atoms < N_atoms){
    //choose every mom atom below CM
    for (int i = 0; i < N_atoms; i++){

      if(pos_mom[i][dir] < pos_CM){
        for (int k = 0; k < nr_atoms; k++){
          flag = 0;
          for (int j = 0; j < 3; j++){
            if (positions[k][j]==pos_mom[i][j])
              flag++;
          }
          if(flag == 3)
            break; 
        }
        if (flag != 3){
          for(int j = 0; j < 3; j++)
            positions[nr_atoms][j] = pos_mom[i][j];
          nr_atoms ++;
          //cout<<"mom below"<<endl;
        }
      }
      if(nr_atoms== N_atoms)
        break;
    }
  }

  if (nr_atoms < N_atoms){
    //choose every dad atom above CM
    for (int i = 0; i < N_atoms; i++){

      if(pos_dad[i][dir] > pos_CM){
        for (int k = 0; k < nr_atoms; k++){
          flag = 0;
          for (int j = 0; j < 3; j++){
            if (positions[k][j]==pos_dad[i][j])
              flag++;
          }
          if(flag == 3)
            break; 
        }
        if (flag != 3){
          for(int j = 0; j < 3; j++)
            positions[nr_atoms][j] = pos_dad[i][j];
          nr_atoms ++;
          //cout<<"dad above"<<endl;
        }
      }
      if(nr_atoms== N_atoms)
        break;
    }
  }

  //f_value = Potencial();

  //delete gRandom;
  /*

  cout<<nr_atoms<<endl;
  cout<<"aaa"<<endl;

  cout << "daughter Atom 1 " << " : " <<posicoes[0][0] << ", " << posicoes[0][1] << ", " << posicoes[0][2] << endl;
  cout << "daughter Atom 2 " << " : " <<posicoes[1][0] << ", " << posicoes[1][1] << ", " << posicoes[1][2] << endl;
  cout << "daughter Atom 3 " << " : " <<posicoes[2][0] << ", " << posicoes[2][1] << ", " << posicoes[2][2] << endl;
  cout << "-------------"<<endl;*/

}
