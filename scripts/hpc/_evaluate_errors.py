#%%
# script_name = str(sys.argv[1])
# job_id = str(sys.argv[2])

#%% import libraries
import glob
import numpy as np


#%% testing "c_download_mrms_mesonet.sh"
# script_name = "c_download_mrms_mesonet.sh"
# job_id = "45124858"
# fpath = "_script_errors/" + script_name

#%% testing "da_cmbn_to_dly_ncs_frmtd_for_RainyDay.sh"
script_name = "da_cmbn_to_dly_ncs_frmtd_for_RainyDay.sh"
job_id = "45250420"
fpath = "_script_errors/" + script_name

#%%
files = glob.glob(fpath + "/" + job_id + "*")

failed_files = []
failed_tasks = []
problem_nodes = []
problems = []
errors = []
files_with_errors = []

for f in files:
    file = open(f, 'r')
    lines = file.readlines()
    for line in lines:
        problem=False
        error=False
        if "error" in line.lower():
            error=True
            problem=True
            tmp_errors = []
            for line in lines:
                tmp_errors.append(line)
            errors.append(tmp_errors)
            files_with_errors.append(f)
        if "DUE TO TIME LIMIT" in line:
            problem=True
            problems.append("FAILED DUE TO TIME LIMIT")
        if "CANCELLED AT" in line:
            problem=True
            problems.append("CANCELLED")
        if problem==True:
            failed_files.append(f)
            details = f.split(job_id)[-1].split("_")
            task = details[1]
            failed_tasks.append(task)
            node = details[-1].split('.')[0]
            node_documented = False
            for n in problem_nodes:
                if n == node:
                    node_documented = True
                    break
            if node_documented == False:
                problem_nodes.append(node)
        if error:
            break

str_tasks = ','.join(str(item) for item in failed_tasks)
str_nodes = ','.join(str(item) for item in problem_nodes)


#%% "c_download_mrms_mesonet.sh" documented excluded nodes so far
excluded_nodes = "udc-aw29-25b,udc-an33-5c0,udc-an33-7c1,udc-aw29-19b,udc-an33-11c1,udc-aw34-3c0,udc-ba26-34c1,udc-aw34-4c0,udc-ba25-32c1,udc-aw29-23a,udc-aw34-19c0,udc-aw34-11c1,udc-aw34-3c1"
excluded_nodes = excluded_nodes.split(',')
num_excluded = len(excluded_nodes)
num_unique = len(np.unique(excluded_nodes))
if num_excluded != num_unique:
    print("Redundant nodes in list")


#%% da_cmbn_to_dly_ncs_frmtd_for_RainyDay inspect errors
with open('_da_errors.txt', 'w') as f:
    for i in np.arange(len(errors)):
        f.write("FILE: {}".format(files_with_errors[i].split(job_id)[1]))
        # print("FILE: {}".format(files_with_errors[i].split(job_id)[1]))
        for l in errors[i]:
            # print(l)
            f.write(l)
        f.write("#################################")


#%% "da_cmbn_to_dly_ncs_frmtd_for_RainyDay" documented excluded nodes so far
tasks_to_rerun = '10,11,12,3,4,5,7,8,9'
excluded_nodes = 'udc-ba26-18,udc-ba27-14,udc-ba26-16,udc-ba26-17'
excluded_nodes = excluded_nodes.split(',')
num_excluded = len(excluded_nodes)
num_unique = len(np.unique(excluded_nodes))
if num_excluded != num_unique:
    print("Redundant nodes in list")
