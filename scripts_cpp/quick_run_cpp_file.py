import os
import sys

###############################################################################
# Specify which file you want to run (w.r.t Mschi_Scwahrz)
runfileRelativeDir = "tests/test_causets_general/"
runfilename = "MMdim_flat_test.cpp"
print(runfileRelativeDir+runfilename)
###############################################################################


cwd = os.getcwd()

# Path to MSci_Schwarzschild_Causets
file_dir = os.path.dirname(os.path.realpath(__file__))
main_dir = os.path.dirname(file_dir)+"/"
runfileDir = main_dir + "scripts_cpp/" + runfileRelativeDir
print(main_dir)

# Directory path of causets_cpp
causets_cpp = main_dir + "scripts_cpp/causets_cpp/"
#causets_cpp = main_dir + "scripts_cpp/causets_cpp_int8_t/"
#print("Using int8_t Cmatrix!")

# Compile commands
compiler_stuff = [
    f"g++ -g '{runfileDir+runfilename}'",
    f"'{causets_cpp}causet.cpp'",
    f"'{causets_cpp}embeddedcauset.cpp'",
    f"'{causets_cpp}sprinkledcauset.cpp'",
    f"'{causets_cpp}shapes.cpp'",
    f"'{causets_cpp}spacetimes.cpp'",
]


# Create Executable
############################################
exec = [ " ",
    "-o",
    #"'output.exe'"
    f"'{runfileDir}{runfilename[:-4]}.exe'"
]

#Optimisation flags
############################################
optimisations = [ " ",
    "-std=c++17",   
    "-Ofast",
    #"-march=native",
    #"-ffast-math",
    "-fopenmp",
    ]

# Include paths
###########################################
from os.path import expanduser
home = expanduser("~")
includes = [" ",
    "-I",
    f"{home}/boost/boost_1_80_0",
    #"/rds/general/user/vh119/home/boost/boost_1_80_0"
    ]

# Compile
############################################
to_render = " ".join(compiler_stuff) + " ".join(exec) + \
            " ".join(optimisations) + " ".join(includes)

print("\nC++ files to compile, flags used, executable produced:")
print(to_render)
#os.system(to_render)
print(".............................Compiled............................!\n")


# Run
############################################
print("     pwd")
os.system("pwd")
print(f'    {runfileDir}{runfilename[:-4]}.exe')
os.system(f'{runfileDir}{runfilename[:-4]}.exe')
print(".............................Run............................!\n")