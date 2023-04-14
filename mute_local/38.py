import subprocess as sub
import pip
import sys
import optuna
import os
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import matplotlib.ticker as mtick
import shutil

saved_trial = True

source_bs = 'best_molecule.bs'
source_mv = 'best_molecule.mv'

dest_bs = '55.bs'
dest_mv = '55.mv'


# define the file name for storing the best overall value
best_overall_file = 'best_overall55.txt'

# check if the file exists and read in the best overall value if it does
if os.path.exists(best_overall_file):
    with open(best_overall_file, 'r') as f:
        best_overall = float(f.read().strip())
else:
	with open(best_overall_file, 'w') as f:
		f.write(str(float('inf')))
	best_overall = float('inf')


def objective(trial):

	fit = 0
	 
	#parameters to optmize and respective intervals
	N_molecules = trial.suggest_int('N_molecules',500,1500)
	survival_rate = trial.suggest_float('survival_rate',0.8,0.99)
	mutation_prob = trial.suggest_float('mutation_prob',0.001,0.1)
	m0 = trial.suggest_float('m0',0.3,0.8)
	L_box = trial.suggest_float('L_box',1.5,2.5)

	#writing the given parameters to .h file
	with open('parameters.h','w') as f:
		f.write('inline unsigned int N_molecules = '+str(N_molecules)+';\n');
		f.write('inline float survival_rate = '+str(survival_rate)+';\n');
		f.write('inline float mutation_prob = '+str(mutation_prob)+';\n');
		f.write('inline float m0 = '+str(m0)+';\n');
		f.write('inline float L_box = '+str(L_box)+';\n');


	#compiling
	process = sub.Popen(["make parallel"],shell =True);
	#wait for the compilation to end
	process.wait();
	# Check if the make command succeeded
	if process.returncode != 0:
		print("Make failed");
		sys.exit(1);

	for i in range(3):#multiple runs for less randomized results
	#running executable and waiting
		process = sub.run('./main',check=True);



		with open('opt.txt', 'r') as f:
			current = float(f.read().strip());

		if current < fit:
			fit = current
			with open(best_overall_file, 'r') as f:
				best_overall = float(f.read().strip())

			if current < best_overall:
				shutil.copy(source_bs,dest_bs)
				shutil.copy(source_mv,dest_mv)
				best_overall = current
				with open(best_overall_file, 'w') as f:
					f.write(str(best_overall))
	return fit



# Create or load the study
study = optuna.create_study()

study.enqueue_trial({'N_molecules': 1000, 'survival_rate': 0.95, 'mutation_prob': 0.01, 'm0': 0.5, 'L_box': 2.})
study.enqueue_trial({
	'N_molecules':968,
	'survival_rate': 0.861766862179341,
	'mutation_prob':0.049580082712412654,
	'm0':0.5831106528245837,
	'L_box':1.7785126964669409})
# Optimize the study
study.optimize(objective, timeout=8*60*60)

trial_df = study.trials_dataframe()
trial_df.to_pickle('my_study.pkl')

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



fig = optuna.visualization.plot_param_importances(study)
fig.show()
fig.write_image("plot_param_importances.pdf")

fig = optuna.visualization.plot_parallel_coordinate(study)
fig.update_layout(yaxis=dict(type='log'))
fig.write_image("plot_parallel_coordinate.pdf")

fig = optuna.visualization.plot_optimization_history(study)
fig.write_image("plot_optimization_history.pdf")

fig = optuna.visualization.plot_slice(study)
fig.write_image("plot_slice.pdf")

fig = optuna.visualization.plot_contour(study)
fig.write_image("plot_contour.pdf")
