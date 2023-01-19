
# Description of Data Folder

The "main" folders with "real" results - with Poiss=True etc. everything working are:

- links
- lambdas
- HRVs

for counting links, lambdas and Hawking Radiation V molecules respectively.

The name of the files is self explanatory. Regarding their content, in each txt file:

- lambdas:

  - Nreps is the number of causet realisation simulated to obtain the statistics.
  - {i}avg and {i}std, for i>0, are the average number and the related standard deviation of lambdas of size {i}.
  - {i}avg and {i}std, for i<0, are the same for: -3 -> radius of outermost element involved in lambda, -2 -> radius of innermost, -1->time coordinates of most in the past.

- HRVs:

  - Nreps is the number of causet realisation simulated to obtain the statistics.
  - {0}avg and {0}std are the average number and the related standard deviation of OPEN HRV molecules.
  - {1}avg and {1}std are same for CLOSE HRV molecules.
  - {i}avg and {i}std, for i<0, are the same for: -3 -> radius of outermost element involved in lambda, -2 -> radius of innermost, -1->time coordinates of most in the past.
