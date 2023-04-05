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
	mutation_prob = trial.suggest_float('mutation_prob',0.001,0.25)
	sex_prob = trial.suggest_float('sex_prob',0.01,0.6)
	mating = trial.suggest_categorical('mating',[0,1])
	m0 = trial.suggest_float('m0',0.05,0.4)
	m_mat = trial.suggest_float('m_mat',0.,1.)

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
study.enqueue_trial({
	'survival_rate':0.9,
	'mutation_prob':0.01,
	'sex_prob':0.5,
	'mating':1,
	'm0':0.3,
	'm_mat':0.5
	})
study.enqueue_trial({
	'survival_rate':0.9,
	'mutation_prob':0.01,
	'sex_prob':0.5,
	'mating':0,
	'm0':0.3,
	'm_mat':0.5
	})
study.enqueue_trial({'survival_rate': 0.8201451367523969, 'mutation_prob': 0.03406594968384061, 'sex_prob': 0.1656358510909437, 'mating': 1, 'm0': 0.24072463591200688, 'm_mat': 0.9432572428301249})
study.enqueue_trial({'survival_rate': 0.8115814339383576, 'mutation_prob': 0.010627441796682041, 'sex_prob': 0.1657359426531127, 'mating': 1, 'm0': 0.23908251017478105, 'm_mat': 0.9457648803181609})
study.enqueue_trial({'survival_rate': 0.8238728622839715, 'mutation_prob': 0.01604364301678423, 'sex_prob': 0.11324159228527274, 'mating': 0, 'm0': 0.253755962568267, 'm_mat': 0.8333342589033529})
study.enqueue_trial({'survival_rate': 0.848584515881648, 'mutation_prob': 0.043098945333826966, 'sex_prob': 0.2958368690535383, 'mating': 0, 'm0': 0.2696843619997009, 'm_mat': 0.9305971949435674})
study.enqueue_trial({'survival_rate': 0.854554816721901, 'mutation_prob': 0.03442794757089047, 'sex_prob': 0.20017020858739346, 'mating': 0, 'm0': 0.2618114198363721, 'm_mat': 0.851492037063979})
study.enqueue_trial({'survival_rate': 0.8592377118832262, 'mutation_prob': 0.09623436350198085, 'sex_prob': 0.21748409772953287, 'mating': 1, 'm0': 0.3542473194078439, 'm_mat': 0.9491510706162831})

study.optimize(objective,n_trials = 600)
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