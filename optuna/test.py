import subprocess as sub
import pip
import sys
import optuna
import os

print (optuna.__version__)

for i in range(5):
	#compiling
	process = sub.Popen(["make parallel"],shell =True)
	#wait for the compilation to end
	process.wait()
	# Check if the make command succeeded
	if process.returncode != 0:
	    print("Make failed")
	    sys.exit(1)

	#running executable
	process = sub.run('./main',check=True)

	print('reading txt file')

	with open('opt.txt', 'r') as f:
	    fit = float(f.read().strip())


	print(fit)

	print('-------------')

	with open('parameters.h','w') as f:
		f.write('inline unsigned int N_molecules = 400;');
		f.write('inline float L_box = 3.;');
		f.write('inline float survival_rate = 0.7;');
		f.write('inline float mutation_prob = 0.15;');
		f.write('inline float sex_prob = 0.20;');

