#include <cstdio>
#include <iostream>
#include <iomanip>
#include "moleculas.h"

using namespace std;

int main()
{
	int N_moleculas = 1;
	int N_atomos = 7;
	int dim_caixa = 1;

	double f_value;

	molecula* mol = new molecula(N_atomos, dim_caixa);

	f_value = mol -> lennard_jones_potential();

	mol -> ~molecula();

	return 0;
}