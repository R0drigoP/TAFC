import subprocess as sub
import pip
import sys
import optuna
import os
import pandas as pd
import numpy as np


def objective(trial):

	calls = 0
	fits = 0
	runs = 0 
	#parameters to optmize and respective intervals
	N_molecules = trial.suggest_int('N_molecules',4,1000)
	survival_rate = trial.suggest_float('survival_rate',0.2,0.97)
	mutation_prob = trial.suggest_float('mutation_prob',0.001,0.25)
	m0 = trial.suggest_float('m0',0.01,0.6)
	L_box = trial.suggest_float('L_box',0.2,2.5)

	#writing the given parameters to .h file
	with open('parameters.h','w') as f:
		f.write('inline unsigned int N_molecules = '+str(N_molecules)+';\n');
		f.write('inline float survival_rate = '+str(survival_rate)+';\n');
		f.write('inline float mutation_prob = '+str(mutation_prob)+';\n');
		f.write('inline float m0 = '+str(m0)+';\n');
		f.write('inline float L_box = '+str(L_box)+';\n');


	#compiling
	process = sub.Popen(["make rodri"],shell =True);
	#wait for the compilation to end
	process.wait();
	# Check if the make command succeeded
	if process.returncode != 0:
		print("Make failed");
		sys.exit(1);

	for i in range(5):#multiple runs for less randomized results
		fit, call = run() #runi and read files
		if fit < -44 :
			calls += call
			fits += fit
			runs +=1

	if runs >0:
		calls = calls/runs
		fits = fits/runs
		return calls

	if runs == 0: #return nan nb of calls so that we dont study the convergence when we haven't achieved the global minimum 
		return 1e6

	else:
		return calls #objective -> minimize calls

def run():
	#running executable and waiting
	process = sub.run('./main',check=True);

	with open('opt.txt', 'r') as f:
		fit = float(f.read().strip());

	with open('calls.txt','r') as f:
		call = int(f.read().strip());

	return fit,call 

study = optuna.create_study()
#study.enqueue_trial({'N_molecules': 318, 'survival_rate': 0.7338862027580697, 'mutation_prob': 0.012303838196609167, 'm0': 0.04945950298427772, 'L_box': 1.1422373438777127})
#study.enqueue_trial({'N_molecules': 63, 'survival_rate': 0.33034961448484323, 'mutation_prob': 0.08778855226404686, 'm0': 0.044911348910886904, 'L_box': 2.1485776666826997})

study.enqueue_trial({'N_molecules': 10, 'survival_rate': 0.8, 'mutation_prob': 0.01, 'm0': 0.02, 'L_box': 0.5})
study.enqueue_trial({'N_molecules': 4, 'survival_rate': 0.9693214336867287, 'mutation_prob': 0.054926374040749036, 'm0': 0.05194502992690725, 'L_box': 0.40354634598545686})
study.enqueue_trial({'N_molecules': 26, 'survival_rate': 0.9678724468035772, 'mutation_prob': 0.018880409901563075, 'm0': 0.1791436839936674, 'L_box': 0.4083221740382169})


study.optimize(objective,timeout=3*60)

best = study.best_params
print(best)

#save best parameters to .h
with open('parameters.h','w') as f:
	f.write('inline unsigned int N_molecules = '+str(best['N_molecules'])+';\n');
	f.write('inline float survival_rate = '+str(best['survival_rate'])+';\n');
	f.write('inline float mutation_prob = '+str(best['mutation_prob'])+';\n');
	f.write('inline float m0 = '+str(best['m0'])+';\n');
	f.write('inline float L_box = '+str(best['L_box'])+';\n');


#save all trials to txt
trials = study.trials_dataframe(attrs = ('number','value','params'))
trials.to_csv('trials.txt',index=False,header=True,sep='\t', float_format='%.7f')

#plots
fig = optuna.visualization.plot_param_importances(study)
fig.show()
fig.write_image("plot_param_importances.pdf")

fig = optuna.visualization.plot_parallel_coordinate(study)
fig.write_image("plot_parallel_coordinate.pdf")

fig = optuna.visualization.plot_optimization_history(study)
fig.write_image("plot_optimization_history.pdf")

fig = optuna.visualization.plot_slice(study)
fig.write_image("plot_slice.pdf")

fig = optuna.visualization.plot_contour(study)
fig.write_image("plot_contour.pdf")

fig = optuna.visualization.plot_contour(study, params=['L_box','m0'])
fig.write_image("plot_contour_L_m0.pdf")

fig = optuna.visualization.plot_contour(study, params=['L_box','mutation_prob'])
fig.write_image("plot_contour_L_mut.pdf")

fig = optuna.visualization.plot_contour(study, params=['L_box','N_molecules'])
fig.write_image("plot_contour_L_mol.pdf")


fig = optuna.visualization.plot_contour(study, params=['L_box','survival_rate'])
fig.write_image("plot_contour_L_surv.pdf")

fig = optuna.visualization.plot_contour(study, params=['N_molecules','m0'])
fig.write_image("plot_contour_mol_m0.pdf")

fig = optuna.visualization.plot_contour(study, params=['N_molecules','mutation_prob'])
fig.write_image("plot_contour_mol_mut.pdf")

fig = optuna.visualization.plot_contour(study, params=['N_molecules','survival_rate'])
fig.write_image("plot_contour_mol_surv.pdf")

fig = optuna.visualization.plot_contour(study, params=['m0','mutation_prob'])
fig.write_image("plot_contour_m0_mut.pdf")

fig = optuna.visualization.plot_contour(study, params=['m0','survival_rate'])
fig.write_image("plot_contour_m0_surv.pdf")

fig = optuna.visualization.plot_contour(study, params=['mutation_prob','survival_rate'])
fig.write_image("plot_contour_mut_surv.pdf")