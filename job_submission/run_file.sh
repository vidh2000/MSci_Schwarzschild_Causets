#!/usr/bin/env bash

echo "=====================================================" 
echo

########################################################################
############################  TO CHANGE ################################
########################################################################

# Files you want to run in 
runfilename="speed_testing_BH.cpp"
#runfilename="test_hello_world.cpp"

# Specify relative path from MSci_Schwarzschild_Causets
# to the .cpp file you want to compile/execute
runfileRelativeDir="scripts_cpp/tests/"

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
python run_cpp_file.py $mainDir$runfileRelativeDir $runfilename

# Go into the directory of the executable file
cd $mainDir$runfileRelativeDir

# Check if inside right directory
#echo "Currently inside -->"
#pwd

# Execute the created .exe program 
#echo ${runfilename::-4}
./${runfilename::-4}".exe"
