# To run the scripts:
- The only relevant hard coded variables and filepaths are in `__directories.sh` and `__utils.py`. Updating these files should be all that's necessary to re-run these scripts.

# Helpful linux commands:
- To combine all the outputs from a single job's worth of outputs or error .out files:
    
    `cd [output or error folder]`
    
    `cat *[job number]*.out > _[job number]_combined.txt`

# Helpful SLURM commands:
- The job ID number is recorded in the output and error filenames

sacct -j [job id #] -o jobid,jobname%20,user,partition,state,start,end,NodeList%60
- this is helpful for getting the cluster’s estimate of when a job will start (start command); you can also submit this for completed jobs to get the wall clock time
- This came in handy when trying to balance memory allocation requests (larger requests would cause a longer wait time) and data chunk size (larger chunks would mean faster computing time)

seff [job id #]
- returns job statistics, e.g., CPU and memory efficiency

# Github notes:
- clone specific branch: https://stackoverflow.com/questions/1911109/how-do-i-clone-a-specific-git-branch/7349740#7349740
    - a simple `git pull` while cd'ed into the repo directory will pull any changes made to that branch
- add a repo within a repo
    - https://stackoverflow.com/questions/1811730/how-do-i-work-with-a-git-repository-within-another-repository