#ifndef __molecule__
#define __molecule__

#include "TRandom3.h"
#include "global.h"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

class molecule{

public:

  //Constructors and Destructor
  molecule(unsigned int n_atoms = 0, float l_box = 0., float mute_prob = 0.);
  molecule(molecule* mom, molecule* dad, double gene_prop, unsigned int n_atoms = 0);
  ~molecule();

  //Fitness
  void Fit();
  
  //Getters
  int Get_Natoms() {return N_atoms;}
  int Get_Dim() {return L_box;}
  double** Get_Pos() {return positions;}
  double Get_Fit() {return fitness;}

  //Setters
  void Set_Pos(double** pos);
  void Set_Fit(double pot);
  
  //Comparator
  static bool LessPot(molecule* mol1, molecule* mol2) {return mol1->Get_Fit() < mol2->Get_Fit();}

  //Mutation and Reproduction
  void Mutate(unsigned int iter, float m0, float alpha, int flag, TRandom3* gRandom = NULL);
  int generate_children3(vector<molecule*> pop, TRandom3* gRandom);
  void Mating(molecule* mom, molecule* dad, double gene_prop = 0.);
  void Mating_Plano3(molecule* mom, molecule* dad, TRandom3* gRandom);

private:
  unsigned int N_atoms;
  float mutation_prob, L_box;
  double **positions;
  double fitness;

};

#endif
