#%% import libraries



#%% hard coding
# associated function: remove_vars
coords_to_delete = ["step", "heightAboveSea", "valid_time"] # do not contain useful information
attrs_to_delete = ['source', 'problems'] # not valid for aggregated timestep

#%% functions
def remove_vars(ds, coords_to_delete, attrs_to_delete):
    '''
    This removes coordinates and attributes that are not useful for the combined MRMS product. These are
    hardcoded in the __utils.py file
    '''

    for crd in coords_to_delete:
        try:
            del ds.coords[crd]
        except: 
            continue
    for att in attrs_to_delete:
        try:
            del ds.attrs[att]
        except:
            continue
    return ds