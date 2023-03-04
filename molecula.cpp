#include "molecula.h"

molecula::molecula(int n_atomos, int dim_caixa){

	N_atomos = n_atomos;

	Dim_caixa = dim_caixa;

	posicoes = new double*[N_atomos];

        for (int i = 0; i < N_atomos; i++) {

            posicoes[i] = new double[3];
        }

    gRandom = new TRandom3(0);

	for (int j = 0; j < N_atomos; j++){
		posicoes[j][0] = gRandom -> Uniform(-Dim_caixa, Dim_caixa); //gera posicao x entre -10 e 10 
		posicoes[j][1] = gRandom -> Uniform(-Dim_caixa, Dim_caixa); //gera posicao y entre -Dim_caixa e Dim_caixa
		posicoes[j][2] = gRandom -> Uniform(-Dim_caixa, Dim_caixa); //gera posicao z entre -10 e 10
		//cout << "Pos: x:"<< posicoes[j][0] << " - y:" << posicoes[j][1] << " - z:"  << posicoes[j][2] << endl;
	}
	
}


molecula::~molecula() {
    //apagar posicoes
    for(int i = 0; i < 3; ++i) { delete[] posicoes[i]; }
    
    delete[] posicoes;
}


//calcula o lennard_jones para cada molecula
double molecula::lennard_jones_potential(){

	float sigma = 1; //(Angstroms)
	float epsilon = 1; //(kJ/mol)

	float d_aux = 0;

	f_value = 0;
	
	//nao sei se esta certo, confirmar com os outros
	for (int j = 0; j < N_atomos; j++){

		for (int k = 0; k < N_atomos; k++){

			if( j < k ){

				d_aux = sqrt( pow((posicoes[j][0] - posicoes[k][0]), 2) + pow((posicoes[j][1] -posicoes[k][1]), 2) + pow((posicoes[j][2] - posicoes[k][2]), 2));

				f_value += 4*epsilon*(pow((sigma/d_aux), 12) - pow((sigma/d_aux), 6));		
			}else{ d_aux = 0;}
		}
	}

	cout << "f_value: "<< f_value << endl;

	return f_value;
}


















