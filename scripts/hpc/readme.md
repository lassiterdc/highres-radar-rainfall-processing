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