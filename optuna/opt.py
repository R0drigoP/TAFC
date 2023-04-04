import subprocess as sub
import pip
import sys
import optuna
import os
import pandas as pd

def objective(trial):
	#N_molecules = trial.suggest_int('N_molecules',400,1000)
	#L_box = trial.suggest_float('L_box',2.5,5.)
	survival_rate = trial.suggest_float('survival_rate',0.5,0.95)
	mutation_prob = trial.suggest_float('mutation_prob',0.001,0.2)
	sex_prob = trial.suggest_float('sex_prob',0.001,0.5)
	mating = trial.suggest_categorical('mating',[0,1])
	m0 = trial.suggest_float('m0',0.05,0.4)
	m_mat = trial.suggest_float('m_mat',0.1,0.9)

	#writing the given parameters to .h file
	with open('parameters.h','w') as f:
		#f.write('inline unsigned int N_molecules = '+str(N_molecules)+';\n');
		#f.write('inline float L_box = '+str(L_box)+';\n');
		f.write('inline float survival_rate = '+str(survival_rate)+';\n');
		f.write('inline float mutation_prob = '+str(mutation_prob)+';\n');
		f.write('inline float sex_prob = '+str(sex_prob)+';\n');
		f.write('inline bool mating = '+str(int(mating))+';\n');
		f.write('inline float m0 = '+str(m0)+';\n');
		f.write('inline float m_mat = '+str(m_mat)+';\n');


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
study.optimize(objective,n_trials = 10)
best = study.best_params
print(best)

with open('parameters.h','w') as f:
	#f.write('inline unsigned int N_molecules = '+str(best['N_molecules'])+';\n');
	#f.write('inline float L_box = '+str(best['L_box'])+';\n');
	f.write('inline float survival_rate = '+str(best['survival_rate'])+';\n');
	f.write('inline float mutation_prob = '+str(best['mutation_prob'])+';\n');
	f.write('inline float sex_prob = '+str(best['sex_prob'])+';\n');
	f.write('inline bool mating = '+str(int(best['mating']))+';\n');
	f.write('inline float m0 = '+str(best['m0'])+';\n');
	f.write('inline float m_mat = '+str(best['m_mat'])+';\n');


trials = study.trials_dataframe(attrs = ('number','value','params'))

trials.to_csv('trials.txt',index=False,header=True,sep='\t', float_format='%.7f')

#print(study.trials)