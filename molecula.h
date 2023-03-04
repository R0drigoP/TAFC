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

	molecula(int n_atomos, int dim_caixa); //
	~molecula(); //

	double lennard_jones_potential();

	private:

	int N_atomos;
	int Dim_caixa;

	double **posicoes;
	double f_value;
	//vector< pair < int, float> > f_value; 
};

#endif