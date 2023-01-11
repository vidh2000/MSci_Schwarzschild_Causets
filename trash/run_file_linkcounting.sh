#!/usr/bin/env bash


# Parameters to be inputted via the link_counting_submitter.sh
mass=$1
N_multiplier=$2
N_reps=$3

########################################################################
############################  TO CHANGE ################################
########################################################################

### Files you want to run in 
runfilename="count_links.cpp"

### Specify relative path from MSci_Schwarzschild_Causets
# to the .cpp file you want to compile/execute
runfileRelativeDir="scripts_cpp/for_submitting_cpp/"


########################################################################
######################## Leave below alone #############################
########################################################################

# Get useful directories
fileDir=$(pwd)"/"
mainDir="$(dirname "$fileDir")/"



#echo $fileDir
#echo $mainDir

# Run the python file which will compile the .cpp file you want
# (the includes,flags... can be adjusted in run_cpp_file.py)
echo Building the executable...
python run_cpp_file.py $mainDir$runfileRelativeDir $runfilename
echo Finished building the executable

# Go into the directory of the executable file
cd $mainDir$runfileRelativeDir


# Execute the created .exe program 
#echo ${runfilename::-4}
./${runfilename::-4}".exe" $mass $N_multiplier $N_reps

