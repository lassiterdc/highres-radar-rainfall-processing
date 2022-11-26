#%%
# script_name = str(sys.argv[1])
# job_id = str(sys.argv[2])

#%% import libraries
import glob
import numpy as np

#%%
def evaluate_errors(script_name, job_id):
    fpath = "_script_errors/" + script_name
    files = glob.glob(fpath + "/" + job_id + "*")

    failed_files = []
    failed_tasks = []
    problem_nodes = []
    problems = []
    errors = []
    files_with_errors = []

    print("inspecting {} files for errors or node failures...".format(len(files)))
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

    num_issues = len(failed_tasks) + len(files_with_errors)

    if num_issues == 0:
        print("There were no errors or node failures in this job!")

    return str_tasks, str_nodes, files_with_errors, errors

def test_for_duplicates(excluded_nodes):
    excluded_nodes = excluded_nodes.split(',')
    num_excluded = len(excluded_nodes)
    num_unique = len(np.unique(excluded_nodes))
    if num_excluded != num_unique:
        print("WARNING: Redundant nodes in list")
    else:
        print("no duplicates.")
    
def write_errors_to_file(script_name, errors):
    f_out = "_errors_script_{}.txt".format(script_name.split('_')[0])
    with open(f_out, 'w') as f:
        for i in np.arange(len(errors)):
            f.write("FILE: {}".format(files_with_errors[i].split(job_id)[1]))
            for l in errors[i]:
                f.write(l)
            f.write("#################################")
#%% "c_download_mrms_mesonet.sh"
script_name = "c_download_mrms_mesonet.sh"
job_id = "45124858"
str_tasks, str_nodes, files_with_errors, errors = evaluate_errors(script_name, job_id)
excluded_nodes = "udc-aw29-25b,udc-an33-5c0,udc-an33-7c1,udc-aw29-19b,udc-an33-11c1,udc-aw34-3c0,udc-ba26-34c1,udc-aw34-4c0,udc-ba25-32c1,udc-aw29-23a,udc-aw34-19c0,udc-aw34-11c1,udc-aw34-3c1"
test_for_duplicates(excluded_nodes)


#%% "da_cmbn_to_dly_ncs_frmtd_for_RainyDay.sh"
script_name = "da_cmbn_to_dly_ncs_frmtd_for_RainyDay.sh"
job_id = "45250420"
str_tasks, str_nodes, files_with_errors, errors = evaluate_errors(script_name, job_id)
excluded_nodes = 'udc-ba26-18,udc-ba27-14,udc-ba26-16,udc-ba26-17'
test_for_duplicates(excluded_nodes)
write_errors_to_file(script_name, errors)

#%% db_resampling_to_hourly_and_daily_timesteps.sh
script_name = "db_resampling_to_hourly_and_daily_timesteps.sh"
job_id = "45275746"
str_tasks, str_nodes, files_with_errors, errors = evaluate_errors(script_name, job_id)
excluded_nodes = ''
# test_for_duplicates(excluded_nodes)
# write_errors_to_file(script_name, errors)


#%% dc_combining_daily_totals_in_annual_netcdfs.sh
script_name = "dc_combining_daily_totals_in_annual_netcdfs.sh"
job_id = "45275750"
str_tasks, str_nodes, files_with_errors, errors = evaluate_errors(script_name, job_id)
excluded_nodes = ''
# test_for_duplicates(excluded_nodes)
# write_errors_to_file(script_name, errors)

#%% ha_generate_annual_statistics_netcdfs.sh
script_name = "ha_generate_annual_statistics_netcdfs.sh"
job_id = "45285799"
str_tasks, str_nodes, files_with_errors, errors = evaluate_errors(script_name, job_id)
#%% hb_generate_annual_statistics_plots.sh
script_name = "hb_generate_annual_statistics_plots.sh"
job_id = "45285800"
str_tasks, str_nodes, files_with_errors, errors = evaluate_errors(script_name, job_id)
#%% i_extract_mrms_data_at_gages.sh
script_name = "i_extract_mrms_data_at_gages.sh"
job_id = "45285801"
str_tasks, str_nodes, files_with_errors, errors = evaluate_errors(script_name, job_id)