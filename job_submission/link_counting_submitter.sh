#!/usr/bin/env bash


###############################################################################
### This file creates necessary files and submits the scripts to the cluster
###############################################################################


# CPP VARIABLES
N_multiplier=1000
N_reps=100

#CLUSTER JOB RESOURCE REQUIREMENTS
ncpus=48
mem=16
runtime="08:00:00" #format: "hh:mm:ss"


# SET MASSES YOU WANT TO SIMULATE
for mass in 1.0 1.5
do 


##############################################################################
##############################################################################
##############       NO NEED TO CHANGE ANYTHING BELOW    #####################
##############################################################################
##############################################################################
##############################################################################


# Directories and other variables
job_submissionsDir="/rds/general/user/vh119/home/MSci_Schwarzschild_Causets/\
                                                            job_submission/"
submitted_jobsDir="${job_submissionsDir}submitted_jobs/"

# Write the run_file.sh which compiles and executes the .cpp script
sh_run_file="run_file_linkcounting.sh"

# Write the submit file based on variables used and run_file.sh to submit
submitfile="${submitted_jobsDir}submit_files/\
       submitme_mass_${mass}_N_multiplier_${N_multiplier}_N_reps_${N_reps}.pbs"

echo "#!/bin/sh" > $submitfile
echo "#PBS -lselect=1:ncpus=${ncpus}:mem=${mem}gb" >> $submitfile
echo "#PBS -lwalltime=${runtime}" >> $submitfile
echo "#PBS -o ${submitted_jobsDir}out_files" >> $submitfile
echo "#PBS -e ${submitted_jobsDir}err_files" >> $submitfile
echo "" >> $submitfile
echo "export OMP_NUM_THREADS=${ncpus}" >> $submitfile
echo "" >> $submitfile
echo "sh ${job_submissionsDir}${sh_run_file} ${mass} ${N_multiplier} ${N_reps}" >> $submitfile
echo "sh ${file}" >> $submitfile
echo "echo Submitted ${file}"

# Submit the job for this mass
echo Submitting $submitfile...
qsub $submitfile

done


