#!/usr/bin/env bash

## DESCRIPTION ## DESCRIPTION ## DESCRIPTION ## DESCRIPTION ## DESCRIPTION ##
###############################################################################
### This file creates necessary files and submits the scripts to the cluster
### IT WAS CREATED BY COPYING THE ONE OF LAMBDAS AND 'LAMBDA'->'HRV'
###############################################################################

# DIRECTORIES
homeDir="${HOME}/MSci_Schwarzschild_Causets/"
job_submissionsDir="${homeDir}job_submission/"
submitted_jobsDir="${job_submissionsDir}submitted_jobs/"
cpp_file_to_run="'count_HRVs.cpp'"  # need 'filename.cpp' inside the string!

# CPP VARIABLES
Rho=5000
N_reps=10


# CLUSTER JOB RESOURCE REQUIREMENTS
ncpus=256
mem=920
runtime="04:00:00" #format: "hh:mm:ss"


# SET MASSES YOU WANT TO SIMULATE
# FIRST ROUND MASSES - 1k ...
# [0.53 0.75 0.92  
# 1.06 1.19 1.30 1.40 1.50 
# 1.59 1.68 1.76 1.84 1.91 1.98 
# 2.05 2.12 2.19 2.25 2.31 2.37 2.43 ] 
# SECOND ROUND MASSES - 1.5k ...
# [0.65 0.84 0.99 1.13 1.24 1.35 1.45 
#  1.55 1.63 1.72 1.88 1.95
#  1.8 
#  2.02 2.09 2.15 2.22 2.28 2.34
#  2.4
#  2.46 ]
counter=0 
# 0.53 0.75 0.92 1.06 1.19 
#1.30 1.40 1.50 1.59 1.68 1.76 1.84 1.91 1.98 
#2.05 2.12 2.19 2.25 2.31 2.37 2.43
#$(seq 2.3 .1 2.5)  2.05 2.12 2.19 2.25 2.31 2.37 - 5 reps 920gb. 1.59 1.68 1.76 1.84 1.91 1.98 - 10 reps 512gb

for mass in 2.37 2.12
do 

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
##############################################################################
#>>>>>>>>>>>>>>>       NO NEED TO CHANGE ANYTHING BELOW    <<<<<<<<<<<<<<<<<<#
##############################################################################
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


#__________________________________________________________________
###### Write the run_file.sh which compiles and executes the .cpp script #####


sh_run_file="${submitted_jobsDir}runfile_files/run_file_HRVcounting_mass_${mass}_N_multiplier_${N_multiplier}_N_reps_${N_reps}.sh"


echo "#!/usr/bin/env bash" > $sh_run_file
echo "" >> $sh_run_file
echo "" >> $sh_run_file
echo "# .cpp file you want to run in" >> $sh_run_file
echo "" >> $sh_run_file 
echo "runfilename=${cpp_file_to_run}" >> $sh_run_file
echo "" >> $sh_run_file
echo "### Specify relative path from MSci_Schwarzschild_Causets" >> $sh_run_file
echo "# to the .cpp file you want to compile/execute" >> $sh_run_file
echo "runfileRelativeDir='scripts_cpp/for_submitting_cpp/'" >> $sh_run_file
echo "" >> $sh_run_file
echo "# Get useful directories" >> $sh_run_file
echo "mainDir='${homeDir}'" >> $sh_run_file
echo "" >> $sh_run_file
echo "# Run the python file which will compile the .cpp file you want" >> $sh_run_file
echo "# (the includes,flags... can be adjusted in run_cpp_file.py)" >> $sh_run_file
echo "echo Building the executable..." >> $sh_run_file
# Copy the run_cpp_file.py into $TMPDIR
echo "python ${job_submissionsDir}run_cpp_file.py" '$mainDir$runfileRelativeDir $runfilename' >> $sh_run_file
echo "echo Finished building the executable" >> $sh_run_file
echo "" >> $sh_run_file
echo "echo Go into the directory of the executable file:" >> $sh_run_file
echo 'cd $mainDir$runfileRelativeDir' >> $sh_run_file
#echo "pwd" >> $sh_run_file
#echo "Files in the directory:"
#echo "ls" >> $sh_run_file
echo "" >> $sh_run_file
echo "" >> $sh_run_file
echo "mass=${mass}" >> $sh_run_file
echo "Rho=${Rho}" >> $sh_run_file
echo "N_reps=${N_reps}" >> $sh_run_file
echo "" >> $sh_run_file
echo "# Execute the copy of the created .exe program" >> $sh_run_file
echo 'cp ${runfilename::-4}".exe" "executables/HRVs_M${mass}_Rho${Rho}_reps${N_reps}.exe"' >> $sh_run_file
#echo './${runfilename::-4}".exe" $mass $N_multiplier $N_reps' >> $sh_run_file
echo '"./executables/HRVs_M${mass}_Rho${Rho}_reps${N_reps}.exe" $mass $Rho $N_reps' >> $sh_run_file
echo "" >> $sh_run_file

###############################################################################
#______________________________________________________________________________
##### Write the submit file based on variables used and run_file.sh to submit
submitfile="${submitted_jobsDir}submit_files/subHRV_mass_${mass}_Rho_${Rho}_N_reps_${N_reps}.pbs"

echo "#!/bin/sh" > $submitfile
echo "#PBS -lselect=1:ncpus=${ncpus}:mem=${mem}gb" >> $submitfile
echo "#PBS -lwalltime=${runtime}" >> $submitfile
echo "#PBS -o ${submitted_jobsDir}out_files" >> $submitfile
echo "#PBS -e ${submitted_jobsDir}err_files" >> $submitfile
echo "" >> $submitfile
echo "export OMP_NUM_THREADS=${ncpus}" >> $submitfile
echo "" >> $submitfile
echo "echo Using ${ncpus} cores" >> $submitfile
echo "echo Using ${mem}GB of RAM" >> $submitfile
echo "" >> $submitfile
echo "echo Running ${sh_run_file}" >> $submitfile
echo "sh ${sh_run_file}" >> $submitfile

#__________________________________________________________________
##### Submit the job for this mass #####
echo Submitting $submitfile
qsub $submitfile

# Increase counter as job was submitted
((counter=counter+1)) 

done

echo "----------------------------------------------------------------------"
echo "Submitted ${counter} jobs for parameters:"
echo "Rho = ${Rho}"
echo "N_reps = ${N_reps}"
echo "ncpus = ${ncpus}"
echo "mem = ${mem}"
echo "runtime = ${runtime}"
echo ""
qstat
echo "========================================================================"


