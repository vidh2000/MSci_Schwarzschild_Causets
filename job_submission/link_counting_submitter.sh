#!/usr/bin/env bash

## DESCRIPTION ## DESCRIPTION ## DESCRIPTION ## DESCRIPTION ## DESCRIPTION ##
###############################################################################
### This file creates necessary files and submits the scripts to the cluster
###############################################################################

# DIRECTORIES
homeDir="/rds/general/user/vh119/home/MSci_Schwarzschild_Causets/"
job_submissionsDir="${homeDir}job_submission/"
submitted_jobsDir="${job_submissionsDir}submitted_jobs/"


# CPP VARIABLES
N_multiplier=400
N_reps=10


# CLUSTER JOB RESOURCE REQUIREMENTS
ncpus=48
mem=16
runtime="00:30:00" #format: "hh:mm:ss"


# SET MASSES YOU WANT TO SIMULATE
counter=0
for mass in 1.0 1.5 2.0
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


sh_run_file="${submitted_jobsDir}runfile_files/run_file_linkcounting_mass_${mass}_N_multiplier_${N_multiplier}_N_reps_${N_reps}.sh"


echo "#!/usr/bin/env bash" > $sh_run_file
echo "" >> $sh_run_file
echo "" >> $sh_run_file
echo "# .cpp file you want to run in" >> $sh_run_file
echo "" >> $sh_run_file 
echo "runfilename='count_links.cpp'" >> $sh_run_file
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
echo "# Go into the directory of the executable file" >> $sh_run_file
echo 'cd $mainDir$runfileRelativeDir' >> $sh_run_file
echo "" >> $sh_run_file
echo "" >> $sh_run_file
echo "# Execute the created .exe program" >> $sh_run_file
echo './${runfilename::-4}".exe" $mass $N_multiplier $N_reps' >> $sh_run_file
echo "" >> $sh_run_file

###############################################################################
#______________________________________________________________________________
##### Write the submit file based on variables used and run_file.sh to submit
submitfile="${submitted_jobsDir}submit_files/submitme_mass_${mass}_N_multiplier_${N_multiplier}_N_reps_${N_reps}.pbs"

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
echo "sh ${sh_run_file}" >> $submitfile
echo "echo Running ${sh_run_file}" >> $submitfile

#__________________________________________________________________
##### Submit the job for this mass #####
echo Submitting $submitfile
qsub $submitfile

# Increase counter as job was submitted
((counter=counter+1)) 

done

echo "----------------------------------------------------------------------"
echo "Submitted ${counter} jobs for parameters:"
echo "N_multiplier = ${N_multiplier}"
echo "N_reps = ${N_reps}"
echo "ncpus = ${ncpus}"
echo "mem = ${mem}"
echo "runtime = ${runtime}"
echo ""
qstat
echo "========================================================================"


