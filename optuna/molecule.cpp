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

void molecule::SelectMutation(unsigned int iter, float m0, float alpha, int flag, TRandom3* gRandom){
  

  double check_if_mute = gRandom->Uniform(0,1);

  if(check_if_mute < mutation_prob){
    nb_of_calls_mute ++;

    double mute_type = gRandom->Uniform(0,1);

    //if (mute_type < 1)
    this -> Mutate(iter, m0, alpha, gRandom);
    //else if(mute_type >1)
      //this -> MutateRotation(gRandom);
    this -> Fit();
  }
  else if (flag == 1)
    this -> Fit();
}


void molecule::MutateRotation(TRandom3* gRandom) {

  this->Fit();

  cout<<"BEFORE "<<this->Get_Fit()<<endl;
  // Select a random rotation angle around each axis
  double theta_x = gRandom->Uniform(-0.1, 0.1);
  double theta_y = gRandom->Uniform(-0.1, 0.1);
  double theta_z = gRandom->Uniform(-0.1, 0.1);

  // Calculate the rotation matrices around each axis
  double R_x[3][3] = {{1, 0, 0}, {0, cos(theta_x), -sin(theta_x)}, {0, sin(theta_x), cos(theta_x)}};
  double R_y[3][3] = {{cos(theta_y), 0, sin(theta_y)}, {0, 1, 0}, {-sin(theta_y), 0, cos(theta_y)}};
  double R_z[3][3] = {{cos(theta_z), -sin(theta_z), 0}, {sin(theta_z), cos(theta_z), 0}, {0, 0, 1}};

  // Calculate the overall rotation matrix
  double R[3][3] = {{0}};
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      for(int k=0; k<3; k++) {
        R[i][j] += R_x[i][k] * R_y[k][j];
      }
    }
  }
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      for(int k=0; k<3; k++) {
        R[i][j] = R[i][j] * R_z[k][j];
      }
    }
  }

  // Apply the rotation to each atom position
  for(int i=0; i<N_atoms; i++) {
    double x = this->Get_Coord(i,0);
    double y = this->Get_Coord(i,1);
    double z = this->Get_Coord(i,2);

    // Store the result of the rotation in a temporary vector
    double temp[3] = {0};
    for(int j=0; j<3; j++) {
      temp[j] = R[j][0] * x + R[j][1] * y + R[j][2] * z;
    }

    // Update the atom position vector with the rotated values
    positions[i][0] = temp[0];
    positions[i][1] = temp[1];
    positions[i][2] = temp[2];
  }  
  this->Fit();

  cout<<"AFTER "<<this->Get_Fit()<<endl;
}

void molecule::Mutate(unsigned int iter, float m0, float alpha, TRandom3* gRandom){

  int atom_to_mutate = (int)gRandom->Uniform(0, N_atoms);
  for(int i=0; i<3; ++i){
    
    double x0 = gRandom->Uniform(-1,1)*m0*L_box;

    double mutation = x0/(1+alpha*iter);
    //cout<<positions[atom_to_mutate][i] <<endl;

    positions[atom_to_mutate][i] += mutation;

    
    //if(positions[atom_to_mutate][i] < 0. || positions[atom_to_mutate][i] > L_box)
    //  positions[atom_to_mutate][i] -= 2*mutation;

  }


}




int molecule::generate_children3(vector<molecule*> pop, TRandom3* gRandom){

  //gRandom = new TRandom3(0);

  double check_if_sex = gRandom->Uniform(0,1);

  //cout<<"chech if sex "<<check_if_sex<<endl;

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
    if(mating_type<m_mat){
      nb_of_calls_mat_plano ++;
      this->Mating_Plano3(pop[mom_index],pop[dad_index], gRandom);
    }

    if(mating_type>m_mat) {
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

void molecule::Mating_Sphere(molecule* mom, molecule* dad, TRandom3* gRandom){

  double **pos_mom = mom -> Get_Pos();
  double **pos_dad = dad -> Get_Pos();

  //min number of atoms of each parent
  const int minAtoms = 2;
  int nAtomsMom = 0;
  int nAtomsDad = 0;
  double r = 0., x=0., y=0., z=0.;

  int flag = 0;
  int nr_atoms = 0;

  //setting positions to 0
  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++)
    positions[i][j] = 0;
  }

  //checking that there are at least 2 atoms for each parent inside the sphere
  while (nAtomsMom <= minAtoms || nAtomsDad <= minAtoms) {
    //randomly selecting the radius 
    r = gRandom->Uniform(0,L_box);

    //center coordinates of sphere using normal dist
    x = gRandom->Gaus();
    y = gRandom->Gaus();
    z = gRandom->Gaus();

    //reset to 0
    nAtomsMom = 0;
    nAtomsDad = 0;

    //loop over atoms
    for (int i = 0; i < N_atoms; i++){
      double dist_mom = sqrt((pos_mom[i][0] - x)*(pos_mom[i][0] - x) + (pos_mom[i][1] - y)*(pos_mom[i][1] - y) + (pos_mom[i][1] - z)*(pos_mom[i][1] - z));
      double dist_dad = sqrt((pos_dad[i][0] - x)*(pos_dad[i][0] - x) + (pos_dad[i][1] - y)*(pos_dad[i][1] - y) + (pos_dad[i][1] - z)*(pos_dad[i][1] - z));
      if (dist_mom <= r) 
        nAtomsMom++;
      if (dist_dad <= r) 
        nAtomsDad++;
    }
  }

  //choose every mom atom inside sphere
  for (int i = 0; i < N_atoms; i++){

    double dist_mom = sqrt((pos_mom[i][0] - x)*(pos_mom[i][0] - x) + (pos_mom[i][1] - y)*(pos_mom[i][1] - y) + (pos_mom[i][1] - z)*(pos_mom[i][1] - z));

    if(dist_mom <= r){
      for(int j = 0; j < 3; j++)
        positions[nr_atoms][j] = pos_mom[i][j];
      nr_atoms ++;
    }
  }

  //choose every dad atom outside sphere
  for (int i = 0; i < N_atoms; i++){
    double dist_dad = sqrt((pos_dad[i][0] - x)*(pos_dad[i][0] - x) + (pos_dad[i][1] - y)*(pos_dad[i][1] - y) + (pos_dad[i][1] - z)*(pos_dad[i][1] - z));
    if(nr_atoms== N_atoms)
      break;
    else if( dist_dad > r){
      for(int j = 0; j < 3; j++)
        positions[nr_atoms][j] = pos_dad[i][j];
      nr_atoms ++;
      //cout<<"dad bellow"<<endl;
    }
  }
  if (nr_atoms < N_atoms){
    //choose every mom outside 
    for (int i = 0; i < N_atoms; i++){
      double dist_mom = sqrt((pos_mom[i][0] - x)*(pos_mom[i][0] - x) + (pos_mom[i][1] - y)*(pos_mom[i][1] - y) + (pos_mom[i][1] - z)*(pos_mom[i][1] - z));

      if( dist_mom > r){
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
      double dist_dad = sqrt((pos_dad[i][0] - x)*(pos_dad[i][0] - x) + (pos_dad[i][1] - y)*(pos_dad[i][1] - y) + (pos_dad[i][1] - z)*(pos_dad[i][1] - z));

      if(dist_dad <= r){
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
}




void molecule::Mating_Plano3(molecule* mom, molecule* dad, TRandom3* gRandom){

  double **pos_mom = mom -> Get_Pos();
  double **pos_dad = dad -> Get_Pos();
  double pos_CM = 0;

  int flag = 0;
  int nr_atoms = 0;

  //gRandom = new TRandom3(0);
  int dir = (int)gRandom -> Uniform(0,3);

  for (int i = 0; i < N_atoms; i++){
    for (int j = 0; j < 3; j++)
    positions[i][j] = 0;
  }

  
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
      for(int j = 0; j < 3; j++)
        positions[nr_atoms][j] = pos_mom[i][j];
      nr_atoms ++;
          //cout<<"mom  above"<<endl;
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
