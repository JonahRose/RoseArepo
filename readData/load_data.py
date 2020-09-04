import h5py
import numpy as np
import os

def get_nfiles(file_name):
    split_file_name = file_name.split("/")
    split_dir_name = split_file_name[:-1]
    dir_name = "/".join(split_dir_name)
    nfiles = len([name for name in os.listdir(dir_name) if split_file_name[-1] in name])
    if nfiles == 0:
        raise NameError(f"Did not find any files in {dir_name} with name {split_file_name[-1]}")
    return nfiles


def get_ndata(file_name, file_type, this=-1, part_key=None):
    if file_type=='group':
        key = 'Ngroups'
    elif file_type=='subgroup':
        key= 'Nsubgroups'
    elif file_type=='part':
        key = 'NumPart'
    else:
        raise NameError("Unknown file type")

    if this == -1: #get total groups rather than a specific file
        with h5py.File(file_name+"0.hdf5",'r') as ofile:
            ngroups = ofile['Header'].attrs[key+"_Total"]
    else:
        with h5py.File(file_name+f"{this}.hdf5",'r') as ofile:
            ngroups = ofile['Header'].attrs[key+"_ThisFile"]

    if part_key != None:
        part_type = int(part_key.split("/")[0][-1])
        return ngroups[part_type]
    return ngroups

#singular key
def get_file_type(key):
    if 'Group' in key:
        return 'group'
    elif 'Subhalo' in key:
        return 'subgroup'
    elif 'PartType' in key:
        return 'part'
    else:
        print(f"Not sure what kind of data you are looking for with {key}")
        return None

#keys should be a list 
#file_name should be an absolute path which end with 
def load_data(path, keys):
    data_dict = dict()
    count_dict = {key:0 for key in keys}
    not_found_dict = {key:[] for key in keys}
    num_files = get_nfiles(file_name)
    for i in range(num_files):
        with h5py.File(f"{file_name}{i}.hdf5", "r") as ofile:
            for key in keys:
                if key not in ofile:
                    not_found_dict[key].append(i)
                    continue
                if i == 0:
                    ar = np.array(ofile[key])
                    file_type = get_file_type(key)
                    if file_type == None:
                        continue
                    ngroups = get_ndata(file_name, this=-1, file_type=file_type, part_key=key)
                    shape = list(ar.shape)
                    shape[0] = ngroups
                    data_dict[key] = np.zeros(shape)
                else:
                    ar = np.array(ofile[key])
                data_dict[key][count_dict[key]:count_dict[key]+ar.shape[0]] = ar
                count_dict[key] += ar.shape[0]
    for key in keys:
        num_not_found = len(not_found_dict[key])
        if num_not_found!=0:
            if num_not_found == num_files:
                print(f"Did not find any data for {key}")
    return data_dict
