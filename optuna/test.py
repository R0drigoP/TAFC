import subprocess as sub
import pip
import sys
import optuna
import os
#print (optuna.__version__)

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



print('end of script file')