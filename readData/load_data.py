import h5py
import numpy as np
import os
import time

def get_nfiles(file_name):
    split_file_name = file_name.split("/")
    split_dir_name = split_file_name[:-1]
    dir_name = "/".join(split_dir_name) + "/"
    nfiles = len([name for name in os.listdir(dir_name) if split_file_name[-1] in name and os.path.isfile(dir_name+name)])
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

    if os.path.isfile(file_name):
        file_name = file_name
    else:
        file_name = f"{file_name}{0 if this == -1 else this}.hdf5"

    if this == -1: #get total groups rather than a specific file
        with h5py.File(file_name,'r') as ofile:
            ngroups = ofile['Header'].attrs[key+"_Total"]
    else:
        with h5py.File(file_name,'r') as ofile:
            ngroups = ofile['Header'].attrs[key+"_ThisFile"]

    if part_key != None:
        part_type = get_part_type(part_key)
        return ngroups[part_type]
    return ngroups

def get_part_type(key):
    return int(key.split("/")[0][-1])

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

#get range of particles for an individual galaxy
def check_contains_gal(path, file_num, gal_file, count, sub_idx, fof_idx):
    if sub_idx + fof_idx == -2:
        return False

    #this will only work if group files are read in before snapshot files
    part_key = None
    if part_key is not None:
        if file_num == gal_file:
            return True
        return False

    elif sub_idx > -1:
        idx = sub_idx
        ngroups = get_ndata(path, this=file_num, file_type='subgroup')
    elif fof_idx > -1:
        idx = fof_idx
        ngroups = get_ndata(path, this=file_num, file_type='group')

    if idx > ngroups + count or idx < count:
        return False
    return True

def get_part_key(file_type, key):
    if file_type=='part':
        part_key = key
    else:
        part_key = None
    return part_key

def get_part_range(part_in_file, data_changed, file_type, offsets, part_type, part_count):

    start = offsets.particles_before_gal[part_type] + data_changed - part_count
    end = start + offsets.particles_in_gal[part_type] - data_changed

    if end > part_in_file: #galaxy is spread accross multiple files 
        end = part_in_file

    data_start = data_changed #start of galaxy particle array
    data_end = data_start + (end-start)

    part_count += end 
    data_changed += data_end - data_start

    return [int(data_start), int(data_end)], [int(start), int(end)], data_changed, part_count


#keys should be a list 
#path should be an absolute path which end with 
#individual galaxy indexing does not support negative indecies
def load_data(path, keys, offsets=None):
    data_dict = dict()
    count_dict = {key:0 for key in keys}
    part_count_dict = {key:0 for key in keys} #counting particles we have encountered before our gal
    data_changed_dict = {key:0 for key in keys} #for acounting for a galaxy spread across many files
    found_dict = {key:0 for key in keys}
    finished_files = {key:-1 for key in keys} #don't double count when loading multiple part types
    num_files = get_nfiles(path)
    done_mass_tab = False
    for i in range(num_files):
        should_continue = False
        file_type = get_file_type(keys[0])
        if offsets is not None and file_type == 'part':
            should_continue = True
            for key in keys:

                part_type = get_part_type(key)
                if i > offsets.part_files[part_type]:
                    continue
                if i < offsets.particles_per_file.shape[0]:
                    num_parts = offsets.particles_per_file[i][part_type]
                    if part_count_dict[key] +num_parts < offsets.particles_before_gal[part_type]:
                        part_count_dict[key] += num_parts
                        finished_files[key]=i
                    else:
                        should_continue = False
                else:
                    should_continue = False

            if should_continue:
                continue

        if os.path.isfile(path):
            file_name = path
        else:
            file_name = f"{path}{i}.hdf5"
        with h5py.File(file_name, "r") as ofile:
            for key in keys:
                if key in ofile:
                    found_dict[key] = 1
                else:
                    continue
                if i == 0 and offsets is None:
                    file_type = get_file_type(key)
                    if file_type == None:
                        continue
                    ar = np.array(ofile[key])
                    part_key = get_part_key(file_type, key)
                    ngroups = get_ndata(path, this=-1, file_type=file_type, part_key=part_key)
                    shape = list(ar.shape)
                    shape[0] = ngroups
                    data_dict[key] = np.zeros(shape)
                           
                if offsets is not None:
                    if file_type == 'part':
                        if i <= finished_files[key]:
                            continue
                        part_type = get_part_type(key)
                        num_parts = ofile["Header"].attrs['NumPart_ThisFile'][part_type]

                        if part_count_dict[key] +num_parts < offsets.particles_before_gal[part_type]:
                            part_count_dict[key] += num_parts
                            continue
                        data_range, ar_range, data_changed_dict[key], part_count_dict[key]  = get_part_range(num_parts, data_changed_dict[key], file_type, offsets, part_type, part_count_dict[key])

                        ar = np.array(ofile[key])
                        if key not in data_dict:
                            shape=list(ar.shape)
                            shape[0] = offsets.particles_in_gal[part_type]
                            data_dict[key] = np.zeros(shape)

                        data_dict[key][data_range[0]:data_range[1]] = ar[ar_range[0]:ar_range[1]]
                    else:
                        if offsets.gal_file == i:
                            ar = np.array(ofile[key])
                            if file_type=='subgroup':
                                data_dict[key] = ar[offsets.sub_idx - offsets.gal_before]
                            elif file_type=='group':
                                data_dict[key] = ar[offsets.fof_idx - offsets.gal_before]
                else:
                    ar = np.array(ofile[key])
                    data_dict[key][count_dict[key]:count_dict[key]+ar.shape[0]] = ar
                    count_dict[key] += ar.shape[0]
    for key in keys:
        if found_dict[key]!=1:
            if key == 'PartType1/Masses':
                if os.path.isfile(path):
                    file_name = path
                else:
                    file_name = f"{path}0.hdf5"
                with h5py.File(file_name, "r") as f:
                    num_parts = f['Header'].attrs['NumPart_Total'][1]
                    mass_tab = np.array(f['Header'].attrs['MassTable'])
                if offsets is None:
                    data_dict[key] = np.zeros(num_parts)+mass_tab[1]
                else:
                    data_dict[key] = np.zeros(int(offsets.particles_in_gal[1])) + mass_tab[1]
            else:
                print(f"Did not find any data for {key}")
    return data_dict
