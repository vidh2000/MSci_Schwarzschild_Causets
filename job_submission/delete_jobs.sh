#!/bin/bash

# Get the list of job IDs
job_ids=$(qstat | awk '{print $1}' | sed '1d')

# Loop through the job IDs and delete each job
for id in $job_ids
do
    qdel "$id"
done

echo "All jobs in the queue have been deleted."