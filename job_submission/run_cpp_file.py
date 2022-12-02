import os

cwd = os.getcwd()

# Path to MSci_Schwarzschild_Causets
file_dir = os.path.dirname(os.path.realpath(__file__))
main_dir = os.path.dirname(file_dir)+"\\"

# Directory path of causets_cpp
causets_cpp = main_dir + "scripts_cpp\\causets_cpp\\"

print(causets_cpp)

# Specify which file you want to run
runfileDir = f"{main_dir}scripts_cpp\\tests\\"
runfilename = "test_hello_world.cpp"
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
exec = [ " ",
    "-o",
    #"'output.exe'"
    f"'{runfileDir}{runfilename[:-4]}.exe'"
]

#Optimisation flags
optimisations = [ " ",
    "-Ofast",
    #"-march=native",
    #"-ffast-math",
    "-fopenmp",
    ]

# Include paths
includes = [" ",
    "-I",
    "'C:\\Program Files\\boost\\boost_1_80_0'",
    #"/rds/general/user/vh119/home/boost/boost_1_80_0",
    ]

# Compile

to_render = " ".join(compiler_stuff) + " ".join(exec) + \
            " ".join(optimisations) + " ".join(includes)

print("To render:")
#print(to_render)
os.system(to_render)

# Run the compiled executable

print("Go into folder:")
#print(runfileDir)
os.system("d:")
os.system(f"cd '{runfileDir}'")

print("To execute:")
#print(f"./{runfilename[:-4]}")
os.system(f"./{exec[1]}")













