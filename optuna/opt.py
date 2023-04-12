import subprocess as sub
import pip
import sys
import optuna
import os
import pandas as pd

def objective(trial):
	N_molecules = trial.suggest_int('N_molecules',20,1000)
	#L_box = trial.suggest_float('L_box',2.5,5.)
	survival_rate = trial.suggest_float('survival_rate',0.2,0.97)
	mutation_prob = trial.suggest_float('mutation_prob',0.001,0.25)
	m0 = trial.suggest_float('m0',0.05,0.6)

	#writing the given parameters to .h file
	with open('parameters.h','w') as f:
		f.write('inline unsigned int N_molecules = '+str(N_molecules)+';\n');
		#f.write('inline float L_box = '+str(L_box)+';\n');
		f.write('inline float survival_rate = '+str(survival_rate)+';\n');
		f.write('inline float mutation_prob = '+str(mutation_prob)+';\n');
		#f.write('inline float sex_prob = '+str(sex_prob)+';\n');
		#f.write('inline bool mating = '+str(int(mating))+';\n');
		f.write('inline float m0 = '+str(m0)+';\n');
		#f.write('inline float m_mat = '+str(m_mat)+';\n');


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

	with open('calls.txt','r') as f:
		calls = int(f.read().strip());

	return fit,calls

study = optuna.create_study(directions=["minimize", "minimize"])
study.enqueue_trial({
	'N_molecules':1000,
	'survival_rate':0.95,
	'mutation_prob':0.01,
	'm0':0.5
	})
study.enqueue_trial({'N_molecules': 715, 'survival_rate': 0.704559811437085, 'mutation_prob': 0.21617117772545139, 'm0': 0.16990717490559384})

study.optimize(objective,n_trials = 25)

best = study.best_params
print(best)

with open('parameters.h','w') as f:
	f.write('inline unsigned int N_molecules = '+str(best['N_molecules'])+';\n');
	#f.write('inline float L_box = '+str(best['L_box'])+';\n');
	f.write('inline float survival_rate = '+str(best['survival_rate'])+';\n');
	f.write('inline float mutation_prob = '+str(best['mutation_prob'])+';\n');
	#f.write('inline float sex_prob = '+str(best['sex_prob'])+';\n');
	#f.write('inline bool mating = '+str(int(best['mating']))+';\n');
	f.write('inline float m0 = '+str(best['m0'])+';\n');
	#f.write('inline float m_mat = '+str(best['m_mat'])+';\n');


trials = study.trials_dataframe(attrs = ('number','value','params'))

trials.to_csv('trials.txt',index=False,header=True,sep='\t', float_format='%.7f')


fig = optuna.visualization.plot_param_importances(study)
fig.show()
fig.write_image("plot_param_importances.pdf")
#print(study.trials)

fig = optuna.visualization.plot_parallel_coordinate(study)
fig.write_image("plot_parallel_coordinate.pdf")

fig = optuna.visualization.plot_optimization_history(study)
fig.write_image("plot_optimization_history.pdf")

fig = optuna.visualization.plot_slice(study)
fig.write_image("plot_slice.pdf")

fig = optuna.visualization.plot_contour(study)
fig.write_image("plot_contour.pdf")