import h5py
import numpy as np
import os

def get_nfiles(file_name):
    split_file_name = file_name.split("/")
    split_dir_name = split_file_name[:-1]
    dir_name = "/".join(split_dir_name)
    nfiles = len([name for name in os.listdir(dir_name) if split_file_name[-1] in name])
    if nfiles == 0:
        print(f"Did not find any files in {dir_name} with name {split_file_name[-1]}")
    return nfiles


def get_ndata(file_name, this=-1, groups=False, subgroups=False, parts=False):
    if not (groups or subgroups or parts):
        raise "Must specify groups, subgroups, or parts"
    
    if groups:
        key = 'Ngroups'
    elif subgroups:
        key= 'Numsubgroups'
    else:
        key = 'NumPart'

    if this == -1: #get total groups rather than a specific file
        with h5py.File(file_name+"0.hdf5",'r') as ofile:
            ngroups = ofile['Header'].attrs[key+"_Total"]
    else:
        with h5py.File(file_name+f"{this}.hdf5",'r') as ofile:
            ngroups = ofile['Header'].attrs[key+"_ThisFile"]       
    return ngroups


def load_FoF(file_name, keys):
    data_dict = dict()
    count_dict = {key:0 for key in keys}
    num_files = get_nfiles(file_name)
    for i in range(num_files):
        with h5py.File(f"{file_name}{i}.hdf5", "r") as ofile:
            for key in keys:
                if key not in ofile:
                    print(f"{key} not found in {file_name+str(i)}.hdf5")
                    continue
                if i == 0:
                    ar = np.array(ofile[key])
                    ngroups = get_ndata(file_name, this=-1, groups=True)
                    shape = list(ar.shape)
                    shape[0] = ngroups
                    data_dict[key] = np.zeros(shape)
                else:
                    ar = np.array(ofile[key])
                data_dict[key][count_dict[key]:count_dict[key]+ar.shape[0]] = ar
                count_dict[key] += ar.shape[0]
    return data_dict

print(load_FoF("/home/j.rose/Projects/my_MW_DAO/IllustrisTNGcode/710DM_center/output/groups_000/fof_tab_000.", ["Group/GroupMass"]))
