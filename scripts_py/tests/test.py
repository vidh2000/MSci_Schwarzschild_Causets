import pandas as pd
from os.path import expanduser
import os


# Home Directory
home = expanduser("~")
# Path
path = os.getcwd()
# specify the input file path
folder = home + f"/MSci_Schwarzschild_Causets/data/test_boundary_vs_density/"

# iterate over the files in the folder
for file in os.listdir(folder):
    if "new" in file:
        continue
    # read the CSV file into a DataFrame
    df = pd.read_csv(os.path.join(folder, file),
                header=None)
    # transpose the data
    df = df.transpose()
    # drop rows that contains NaN data
    df = df.dropna()
    # write the transposed data to the original file
    nn = file+"_new"
    data = df.to_string(index=False, header=False, justify='left',
            max_colwidth=20)
    # write the data to the original file
    with open(os.path.join(folder, file), 'w') as f:
        f.write(data)