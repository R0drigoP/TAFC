#ifndef __Molecule__
#define __Molecule__

class Molecule{
private:
  unsigned short int N_atoms, L_box;
  float mutation_rate;
  double** positions;
  double fit_value;
  
public:
  //Constructors and Destructor
  Molecule(unsigned short int n_mol = 0, unsigned short int l_box = 0, float mute_rate = 0.);
  ~Molecule();

  //Getters
  unsigned short int Get_Natoms() {return N_atoms;}
  unsigned short int Get_Lbox() {return L_box;}
  double** Get_Pos() {return positions;}
  double Get_Fit() {return fit_value;}
  
  //Setters
  void Set_Pos(double** new_pos);
  
  //Operations
  static bool LessPot(Molecule* mol1, Molecule* mol2) {return mol1->Get_Fit() < mol2->Get_Fit();}
  //static bool LessPot(Molecule& mol1, Molecule& mol2) {return mol1.Get_Fit() < mol2.Get_Fit();}
  
  //Calculations
  void Fit();
  void Mutate();
};

#endif
