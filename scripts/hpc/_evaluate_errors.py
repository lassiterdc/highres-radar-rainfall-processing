#%%
import glob
script_name = str(sys.argv[1])
job_id = str(sys.argv[2])

#%% testing
script_name = "c_download_mrms_mesonet.sh"
fpath = "_script_errors/" + script_name
job_id = "44848142"
#%%
files = glob.glob(fpath + "/" + job_id + "*")

failed_files = []
failed_tasks = []
problem_nodes = []

for f in files:
    file = open(f, 'r')
    lines = file.readlines()
    for line in lines:
        if "DUE TO TIME LIMIT" in line:
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

str_tasks = ','.join(str(item) for item in failed_tasks)
str_nodes = ','.join(str(item) for item in problem_nodes)