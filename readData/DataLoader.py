import readData.load_data as load
from readData.Offsets import Offsets

import h5py
import numpy as np
import os
import time

class DataLoader():

    def __init__(self, path, snap_num, part_types=-1, keys=[], sub_idx=-1, fof_idx=-1, sub_offsets=None, fof_offsets=None):
        self._check_input(path, part_types, keys) #TODO:currently implemented in parts of many functions

        self.part_types = self._fix_part_types(part_types)
        self.snap_num = str(snap_num).zfill(3)
        
        self.sub_idx = sub_idx
        self.fof_idx = fof_idx

        self.path = path
        self.snap_path = ''
        self.group_path = ''
        self.get_paths()

        self.boxsize = None
        self.time = None
        self.redshift = None
        self.h = None
        self.masstable = None
        self.num_parts = None
        self.num_subhalos = None
        self.num_halos = None
        self.num_part_files = None
        self.num_grp_files = None
        self._get_header_info()

        #create Offset objects if needed (only want a specific galaxy/group)
        self.part_offsets = None
        self.sub_offsets = None
        self.fof_offsets = None
        self._assign_offsets()

        #Change 'GroupMass' -> 'Groups/GroupMass'
        self.keys = [] 
        if type(keys)==type([]):
            if len(keys) > 0:
                self.keys = self.get_correct_keys(keys)
        elif type(keys) == type('string'):
                self.keys = self.get_correct_keys([keys])

        self.file_keys = self.get_file_key_pairs()

        self.data = dict()
        self._load_data()

    def __repr__(self):
        return_str = ''
        return_str += f"Path: {self.path}\n"
        return_str += f"Snap Number: {self.snap_num}\n"
        return_str += f"Part Types: {self.part_types}\n"
        return_str += f"Keys: {self.data.keys()}"
        return return_str

    def __contains__(self, attr):
        key = self.get_correct_keys([attr])[0]
        return key in self.data

    def __getitem__(self, attr):
        key = self.get_correct_keys([attr])[0]
        if key not in self.data:
            raise KeyError(f"Did not load {attr}")
        return self.data[key]

    def _assign_offsets(self):

        if self.fof_idx ==-1:
            if self.sub_idx != -1: #only sub
                self.sub_offsets = Offsets(self.group_path, self.snap_path, self.sub_idx, -1, self.num_grp_files, self.num_part_files)   
        else: #has fof
            self.fof_offsets = Offsets(self.group_path, self.snap_path, -1, self.fof_idx, self.num_grp_files, self.num_part_files) 
            if self.sub_idx != -1: #has both
                self.sub_offsets = Offsets(self.group_path, self.snap_path, self.sub_idx, self.fof_idx, self.num_grp_files, self.num_part_files)

        #figure out which Offset particles should use (default to subfind)
        if self.sub_idx == -1:
            self.part_offsets = self.fof_offsets
        else:
            self.part_offsets = self.sub_offsets

        return


    def _load_data(self):
        self.data = dict()
        for typ in self.file_keys.keys():
            if typ == 'subhalo':
                data_typ = load.load_data(self.group_path, self.file_keys[typ], self.sub_offsets)
            elif typ == 'group':
                data_typ = load.load_data(self.group_path, self.file_keys[typ], self.fof_offsets)
            elif typ == 'part':
                data_typ = load.load_data(self.snap_path, self.file_keys[typ], self.part_offsets)
            else:
                raise NameError("Need to propagate changes to self.file_keys") 
            for key in data_typ:
                self.data[key] = data_typ[key]
        return

    #take the path input and create snap and, group paths
    #snap and group paths end so just the file number is required
    #eg. snapdir_###/snap_###.
    def get_paths(self):
        if 'hdf5' in self.path:
            self.snap_path = self.path
            self.group_path = self.path
            return
        if self.path[-1] != "/":
            self.path += "/"
        if 'output' in os.listdir(self.path):
            self.path += 'output/'
        indir = os.listdir(self.path)

        snap_dirs = [name for name in indir if 'snap' in name] # and not os.path.isfile(self.path + name)]
        if len(snap_dirs) > 0:
            path_list = [name for name in snap_dirs if self.snap_num in name]
            if len(path_list)==0:
                raise ValueError(f"Snap {self.snap_num} does not appear to be present in {self.path}")
            select_path = path_list[0]
            if len(path_list)>1: #check if we find more than one possible snapshot file
                for path in path_list:
                    if 'snapdir' in path or '.hdf5' in path:
                        select_path = path
                print(f"Found more than one possible snapshot folder, I am choosing to use: {select_path}")
            self.snap_path = self.path + select_path
            if not os.path.isfile(self.snap_path): #check if the snapshot contains more than one file
                self.snap_path += '/'
                snap_files = os.listdir(self.snap_path)
                self.snap_path += [name for name in snap_files if '.hdf5' in name][0].split('.')[0] +'.'

        #do the same for the group files
        group_dirs = [name for name in indir if 'group' in name or 'fof' in name]
        if len(group_dirs) > 0:
            group_list = [name for name in group_dirs if self.snap_num in name]
            if len(group_list) == 0:
                raise ValueError(f'Group {self.snap_num} does not appear to be present in {self.path}')
            select_group = group_list[0]
            if len(group_list) > 1:
                for path in group_list:
                    if 'groups' in path or '.hdf5' in path:
                        select_group = path 
            self.group_path = self.path + select_group 
            if not os.path.isfile(self.group_path):
                self.group_path += '/'
                group_files = os.listdir(self.group_path)
                self.group_path += [name for name in group_files if '.hdf5' in name][0].split('.')[0] + '.'
            else:
                group_files = [self.group_path]

            if 'subhalo' not in group_files[0] and self.sub_idx != -1:
                raise NameError("Trying to get a subhalo, but no subfind data is present")

        if len(group_dirs) == 0 and len(snap_dirs) == 0:
            raise NameError(f"No snapshot data found at {self.path}")

        return

    def _fix_part_types(self, part_types):
        ty = type(part_types)
        if ty == type([]):
            return part_types
        elif ty == type(1) or ty == type(1.0):
            if part_types == -1:
                return [0,1,2,3,4,5]
            if part_types >= 5 or part_types < 0:
                print(f"Did you mean PartType{part_types}?")
            return [part_types]
        else:
            raise NameError("PartType not understood")

    def _get_header_info(self):
        print(self.snap_path)
        if os.path.isfile(self.snap_path):
            file_name = self.snap_path
        else:
            file_name = self.snap_path + '0.hdf5'
        with h5py.File(file_name, "r") as ofile:
            pheader = ofile['Header']
            self.boxsize = float(pheader.attrs['BoxSize'])
            self.num_parts = pheader.attrs['NumPart_Total']
            self.time = pheader.attrs['Time']
            self.redshift = pheader.attrs['Redshift']
            self.num_part_files = int(pheader.attrs['NumFilesPerSnapshot'])
            self.h = float(pheader.attrs['HubbleParam'])
            self.masstable = np.array(pheader.attrs['MassTable'])

        if os.path.isfile(self.group_path):
            file_name = self.group_path
        else:
            file_name = self.group_path + '0.hdf5'
        with h5py.File(file_name, "r") as ofile:
            gheader = ofile['Header']
            
            if 'Nsubgroups_Total' in gheader.attrs:
                self.num_subhalos = int(gheader.attrs['Nsubgroups_Total'])
                self.num_halos = int(gheader.attrs['Nsubgroups_Total'])
                self.num_grp_files = int(gheader.attrs['NumFiles'])

                if self.sub_idx > self.num_subhalos:
                    raise KeyError(f"This snapshot only has {self.num_subhalos} subhalos, tried to access subhalo {self.sub_idx}")
                if self.fof_idx > self.num_halos:
                    raise KeyError(f"This snapshot only has {self.num_halos} fof groups, tried to access fof group {self.fof_idx}")

        return 

    def get_file_key_pairs(self):
        file_keys = dict() # {'group':[], 'subhalo':[], 'part':[]}
        for key in self.keys:
            if 'Group' in key:
                file_keys.setdefault('group', []).append(key)
            elif 'Subhalo' in key:
                file_keys.setdefault('subhalo', []).append(key)
            else:
                file_keys.setdefault('part', []).append(key)
        return file_keys

    def get_correct_keys(self, input_keys):
        corrected_keys = []
        for key in input_keys:
            if len(key.split('/')) == 2: #if key is already correct
                corrected_keys.append(key)
                continue
            if 'Group' in key:
                corrected_keys.append(f'Group/{key}')
            elif 'Subhalo' in key:
                corrected_keys.append(f'Subhalo/{key}')
            else:
                for i in self.part_types:
                    corrected_keys.append(f"PartType{i}/{key}")
        return corrected_keys

    def get_halo(self, part_types=-1, sub_idx=-1, fof_idx=-1, keys=[]):

        new_halo = DataLoader(self.path, self.snap_num, part_types, sub_idx=sub_idx, fof_idx=fof_idx, sub_offsets=self.sub_offsets, fof_offsets=self.fof_offsets)

        return new_halo

    #TODO
    def _check_input(self, path, part_types, keys):
        return
