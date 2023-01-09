#!/usr/bin/env bash

# Variables for the cpp script

# CPP VARIABLES
N_multiplier=1000
N_reps=100

#CLUSTER JOB RESOURCE REQUIREMENTS
ncpus=48
mem=16



# SET MASSES YOU WANT TO SIMULATE
for mass in 1.0 1.5
do 


##############################################################################
##############################################################################
##############          NO NEED TO CHANGE ANYTHING BELOW #####################
##############################################################################
##############################################################################
##############################################################################

submitfile="submit_files/submitme_mass_${mass}_N_multiplier_${N_multiplier}_N_reps_${N_reps}.pbs"
# Write the submit file for the specific mass
echo "#!/bin/sh" > $submitfile
echo "#PBS -lselect=1:ncpus=${ncpus}:mem=${mem}gb" >> $submitfile
echo "#PBS -lwalltime=08:00:00" >> $submitfile
echo "#PBS -o /rds/general/user/vh119/home/MSci_Schwarzschild_Causets/job_submission/test_outfolder" >> $submitfile
echo "#PBS -e /rds/general/user/vh119/home/MSci_Schwarzschild_Causets/job_submission/test_errfolder" >> $submitfile
echo "" >> $submitfile
echo "export OMP_NUM_THREADS=${ncpus}" >> $submitfile
echo "" >> $submitfile
echo "file=/rds/general/user/vh119/home/MSci_Schwarzschild_Causets/job_submission/run_file_linkcounting.sh ${mass} ${N_multiplier} ${N_reps}" >> $submitfile
echo "sh ${file}" >> $submitfile

# Submit the job for this mass
echo Submitting $submitfile...
qsub $submitfile

done


