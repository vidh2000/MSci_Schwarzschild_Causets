# Instructions For Submission

From inside job_submission folder, `sh \{name\}_submitter.sh` submits files for analysis, with name being either

- linkcounting
- lambda_counting
- HRVs_counting

In general, the `{name\}_submitter.sh` file has a variable `cpp_file_to_run`: this is the file that needs to be run. The `{name\}_submitter.sh` file creates the sh file running `cpp_file_to_run`. 

Then stuff will appear (thanks to the other files created).

Other things are "old", before the proper job submission scripts were written
which submit 1 job per mass.

## Requirements

Requires miniconda and boost. The latter works with version 1.80.0, available as <https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz/download>. No idea for the other versions. These have to be installed as:

- homedir/miniconda
- homedir/boost/boost_1_80_0

It is also required to configure github log. This must be done by commands

- `git config --global user.email "myeamail@whatever.smth"`
- `git config --global user.name "myusername"`

## Results

The results are then written into the text file with avg and std of the links/lambdas/HRV molecules
counted for specific causet generation configuration.
