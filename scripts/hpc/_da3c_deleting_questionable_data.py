#%% import libraries
from glob import glob
import os
from __utils import *

#%% work
fldr_out_nc = "/scratch/dcl3nd/highres-radar-rainfall-processing/data/mrms_nc_preciprate_fullres_dailyfiles_constant_tstep/"
#%% end work
fldr_out_nc = str(sys.argv[1])

lst_f_nc = glob(fldr_out_nc + "*.nc")
lst_f_nc.sort()

for f in lst_f_nc:
    if "qaqc" in f: # skip qaqc netcdfs
        continue
    # determine date associated with the netcdf file
    str_date = f.split("/")[-1].split(".")[0]
    date = pd.to_datetime(str_date)
    # if that date is in the list of dates to omit for rainyday by renaming it
    if date in dates_to_omit_for_rainyday:
        fnew = f.split(str_date)[0] + "omit_from_rainyday_based_on_qaqc_" + str_date + ".netcdf"
        print("renaming {} to {}".format(f, fnew))
        os.rename(f, fnew)