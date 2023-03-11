#include <cstdlib>
#include <time.h>
#include <cmath>

#include "TRandom3.h"

#include "Molecule.h"

Molecule::Molecule(unsigned short int n_atoms, unsigned short int l_box, float mute_rate) : N_atoms(n_atoms), L_box(l_box), mutation_rate(mute_rate), fit_value(0.){
  positions = new double*[N_atoms];

  for(int i=0; i<N_atoms; ++i)
    positions[i] = new double[3];

  gRandom = new TRandom3(0);

  for(int i=0; i<n_atoms; ++i){
    for(int j=0; j<3; ++j){
      double x = gRandom->Uniform(-l_box/2, l_box/2); //random number between -L_box/2 and L_box/2
      positions[i][j] = x;
    }
  }
}

Molecule::~Molecule(){
  for(int i=0; i<N_atoms; ++i)
    delete [] positions[i];
  delete [] positions;
}

//Setters
void Molecule::Set_Pos(double** new_pos){
  for(int i=0; i<N_atoms; ++i){
    for(int j=0; j<3; ++j)
      positions[i][j] = new_pos[i][j];
  }
}

//Calculations
void Molecule::Fit(){
  double f = 0.;

  for(int i=0; i<N_atoms; ++i){
    for(int j=i+1; j<N_atoms; ++j){
      double r = 0.;
      for(int k=0; k<3; ++k)
	r += (positions[i][k] - positions[j][k])*(positions[i][k] - positions[j][k]);
      r = sqrt(r);

      f += pow(1./r, 12) - pow(1./r, 6);
    }
  }

  fit_value = f;
}


void Molecule::Mutate(){
  gRandom = new TRandom3(0);
  double check_if_mute = gRandom->Uniform(0,1);

  if(check_if_mute < mutation_rate){
    int atom_to_mutate = (int)gRandom->Uniform(0, N_atoms);
    for(int i=0; i<3; ++i){
      double mutation = gRandom->Uniform(-1,1)*0.01*L_box;
      positions[atom_to_mutate][i] += mutation;
      if(positions[atom_to_mutate][i] > L_box)
	positions[atom_to_mutate][i] -= 2*mutation;
    }
  }
}

