#!/usr/bin/env bash


########################################################################
############################  TO CHANGE ################################
########################################################################

### Files you want to run in 
#runfilename="link_counting_4D_futmatrix.cpp"
runfilename="generate_causet_lambdas.cpp"
#runfilename="speed_testing_BH.cpp"
#runfilename="test_hello_world.cpp"

### Specify relative path from MSci_Schwarzschild_Causets
# to the .cpp file you want to compile/execute
runfileRelativeDir="scripts_cpp/for_submitting_cpp/"




########################################################################
######################## Leave below alone #############################
########################################################################

# Get useful directories
homeDir="${HOME}/MSci_Schwarzschild_Causets/"
job_submissionsDir="${homeDir}job_submission/"
submitted_jobsDir="${job_submissionsDir}submitted_jobs/"

mainDir=$homeDir

# Run the python file which will compile the .cpp file you want
# (the includes,flags... can be adjusted in run_cpp_file.py)
echo Building the executable...
python ${job_submissionsDir}run_cpp_file.py $mainDir$runfileRelativeDir $runfilename
echo Finished building the executable

# Go into the directory of the executable file
cd $mainDir$runfileRelativeDir


# Execute the copy of the created .exe program
cp ${runfilename::-4}".exe" "executables/${runfilename::-4}.exe"
./executables/${runfilename::-4}".exe" 
