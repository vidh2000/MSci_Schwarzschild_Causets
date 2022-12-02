import os

cwd = os.getcwd()

# Path to MSci_Schwarzschild_Causets
file_dir = os.path.dirname(os.path.realpath(__file__))
main_dir = os.path.dirname(file_dir)
print("Main dir:\n", main_dir)

# Directory path of causets_cpp
causets_cpp = main_dir + "/scripts_cpp/causets_cpp"

# Specify which file you want to run
runfile_name = main_dir #+ ..........

# Compile commands
compiler_stuff = [
    f"g++ -g {runfile_name}",
    f"{causets_cpp}/causet.cpp",
    f"{causets_cpp}/embeddedcauset.cpp",
    f"{causets_cpp}/sprinkledcauset.cpp",
    f"{causets_cpp}/shapes.cpp",
    f"{causets_cpp}/spacetimes.cpp",
]

# Create Executable
exec = [
    "-o",
    f"{runfile_name[:-4]}"
]

#Optimisation flags
optimisations = [
    "-Ofast",
    #"-march=native",
    #"-ffast-math",
    "-fopenmp",
    ]

# Include paths
includes = [
    "-I",
    "/rds/general/user/vh119/home/boost/boost_1_80_0",
    ]

# Compile

to_render = " ".join(compiler_stuff) + " ".join(exec) + \
            " ".join(optimisations) + " ".join(includes)
os.system(to_render)

# Run the compiled executable
os.system(f"./{exec[1]}")













