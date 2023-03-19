#ifndef __molecula__
#define __molecula__

#include <cstdio>
#include <iostream>
#include <iomanip>
#include "TRandom3.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

class molecula{
public:

  //Constructors and Destructor
  molecula(int n_atomos=0, int dim_caixa=0, double mutation_rate=0.);
  molecula(const molecula&);
  molecula(molecula* mom, molecula* dad, double gene_prop, int n_atomos=0);
  ~molecula();

  //Getters
  int Get_Natomos() {return N_atomos;}
  int Get_Dim() {return Dim_caixa;}
  double** Get_Pos() {return posicoes;}
  double Get_Pot() {return f_value;}

  //Setters
  void Set_Pos(double** pos) {posicoes = pos;}
  
  //Comparator
  static bool LessPot(molecula* mol1, molecula* mol2) {return mol1->Get_Pot() < mol2->Get_Pot();}
  void operator=(const molecula&); //copia os vetores de uma para outra

  //Calculations
  double Potencial();
  double OtherPotential();

  void Mating(molecula* mom, molecula* dad, double gene_prop=0);
  void Mating_Plano(molecula* mom, molecula* dad, double gene_prop=0);
  void Mutate();

  void Mutate_1Atom();
private:
  int N_atomos;
  int Dim_caixa;

  double **posicoes;
  double f_value;
  double mute_rate;
  //vector< pair < int, float> > f_value; 
};

#endif