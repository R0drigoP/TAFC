import subprocess as sub
import pip
import sys
import optuna
import os

def objective(trial):
	N_molecules = trial.suggest_int('N_molecules',30,1000)
	L_box = trial.suggest_float('L_box',2.,10.)
	survival_rate = trial.suggest_float('survival_rate',0.1,0.9)
	mutation_prob = trial.suggest_float('mutation_prob',0.001,0.5)
	sex_prob = trial.suggest_float('sex_prob',0.,0.5)
	mating = trial.suggest_categorical('mating',[0,1])

	#writing the given parameters to .h file
	with open('parameters.h','w') as f:
		f.write('inline unsigned int N_molecules = '+str(N_molecules)+';\n');
		f.write('inline float L_box = '+str(L_box)+';\n');
		f.write('inline float survival_rate = '+str(survival_rate)+';\n');
		f.write('inline float mutation_prob = '+str(mutation_prob)+';\n');
		f.write('inline float sex_prob = '+str(sex_prob)+';\n');
		f.write('inline bool mating = '+str(int(mating))+';\n');

	#compiling
	process = sub.Popen(["make parallel"],shell =True);
	#wait for the compilation to end
	process.wait();
	# Check if the make command succeeded
	if process.returncode != 0:
		print("Make failed");
		sys.exit(1);

	#running executable and waiting
	process = sub.run('./main',check=True);

	with open('opt.txt', 'r') as f:
		fit = float(f.read().strip());

	return fit

study = optuna.create_study()
study.optimize(objective,n_trials=15)
print(study.best_params)