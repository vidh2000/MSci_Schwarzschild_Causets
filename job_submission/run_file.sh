#!/usr/bin/env bash

echo "=====================================================" 
echo


fileDir=$(pwd)"/"
mainDir="$(dirname "$fileDir")/"

# Specify relative path from MSci_Schwarzschild_Causets
# to the .cpp file you want to compile/execute
runfilename="speed_testing_BH.cpp"
runfileRelativeDir="scripts_cpp/tests/"

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