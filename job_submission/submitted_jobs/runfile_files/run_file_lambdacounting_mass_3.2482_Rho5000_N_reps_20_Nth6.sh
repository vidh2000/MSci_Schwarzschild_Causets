#!/usr/bin/env bash


# .cpp file you want to run in

runfilename='count_lambdas_4D.cpp'

### Specify relative path from MSci_Schwarzschild_Causets
# to the .cpp file you want to compile/execute
runfileRelativeDir='scripts_cpp/for_submitting_cpp/'

# Get useful directories
mainDir='/rds/general/user/vh119/home/MSci_Schwarzschild_Causets/'

# Run the python file which will compile the .cpp file you want
# (the includes,flags... can be adjusted in run_cpp_file.py)
echo Building the executable...
python /rds/general/user/vh119/home/MSci_Schwarzschild_Causets/job_submission/run_cpp_file.py $mainDir$runfileRelativeDir $runfilename
echo Finished building the executable

echo Go into the directory of the executable file:
cd $mainDir$runfileRelativeDir


mass=3.2482
Rho=5000
N_reps=20
Nth=6

# Execute the copy of the created .exe program
cp ${runfilename::-4}".exe" "executables/lambdas_M${mass}_Rho${Rho}_reps${N_reps}_Nth${Nth}.exe"
"./executables/lambdas_M${mass}_Rho${Rho}_reps${N_reps}_Nth${Nth}.exe" $mass $Rho $N_reps

