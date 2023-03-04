#include <cstdio>
#include <iostream>
#include <iomanip>
#include "moleculas.h"

using namespace std;

int main()
{
	int N_moleculas = 4;
	int N_atomos = 2;
	int N_dimensoes = 3;

	moleculas* mol = new moleculas(N_moleculas, N_atomos, N_dimensoes);

	cout << "ola" << endl;
	mol -> lennard_jones_function();
	mol -> melhores_moleculas();

	mol -> ~moleculas();

	return 0;
}