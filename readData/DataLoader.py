import readData.load_data as load
import numpy as np
import os


class DataLoader():

    def __init__(self, path, snap_num, part_types, keys):
        self._check_input(path, part_types, keys)

        self.part_types = self._fix_part_types(part_types)
        self.snap_num = str(snap_num).zfill(3)
        
        self.path = path
        self.snap_path = ''
        self.group_path = ''
        self.get_paths()
        
        self.pt1_mass = self.get_pt1_mass()

        #Change 'GroupMass' -> 'Groups/GroupMass'
        self.keys = [] #TODO check what this does 
        if type(keys)==type([]):
            if len(keys) > 0:
                self.keys = self.get_correct_keys(keys)
        elif type(keys) == type('string'):
                self.keys = self.get_correct_keys([keys])

        self.file_keys = self.get_file_key_pairs()

        self.data = dict()
        self._load_data()

    #TODO
    def __repr__(self):
        return ''

    #TODO: check if key is in data
    def __getitem__(self, attr):
        key = self.get_correct_keys([attr])
        return self.data[key[0]]

    def _load_data(self):
        self.data = dict()
        for typ in self.file_keys.keys():
            if typ in ['group', 'subhalo']:
                data_typ = load.load_data(self.group_path, self.file_keys[typ])
            elif typ == 'part':
                data_typ = load.load_data(self.snap_path, self.file_keys[typ])
            else:
                raise NameError("Need to propagate changes to self.file_keys") 
            for key in data_typ:
                self.data[key] = data_typ[key]
        return

    #take the path input and create snap and, group paths
    #snap and group paths end so just the file number is required
    #eg. snapdir_###/snap_###.
    def get_paths(self):
        if self.path[-1] != "/":
            self.path += "/"
        if 'output' in os.listdir(self.path):
            self.path += 'output/'
        indir = os.listdir(self.path)

        snap_dirs = [name for name in indir if 'snap' in name and not os.path.isfile(self.path + name)]
        if len(snap_dirs) > 0:
            self.snap_path = self.path + [name for name in snap_dirs if self.snap_num in name][0] + '/'
        snap_files = os.listdir(self.snap_path)
        self.snap_path += [name for name in snap_files if '.hdf5' in name][0].split('.')[0] +'.'

        group_dirs = [name for name in indir if 'group' in name and not os.path.isfile(self.path + name)]
        if len(group_dirs) > 0:
            self.group_path = self.path + [name for name in group_dirs if self.snap_num in name][0] + '/'
        group_files = os.listdir(self.group_path)
        self.group_path += [name for name in group_files if '.hdf5' in name][0].split('.')[0] + '.'

        return

    def _fix_part_types(self, part_types):
        ty = type(part_types)
        if ty == type([]):
            return part_types
        elif ty == type(1) or ty == type(1.0):
            if part_types >= 5 or part_types < 0:
                print(f"Did you mean PartType{part_types}?")
            return [part_types]
        else:
            raise NameError("PartType not understood")

    def get_pt1_mass(self):
        if 1 in self.part_types:
            return 0
        with h5py.File(self.snap_path + '0.hdf5') as ofile:
            masses = ofile['Header'].attrs['MassTable']
        return masses[1]

    def get_file_key_pairs(self):
        file_keys = {'group':[], 'subhalo':[], 'part':[]}
        for key in self.keys:
            if 'Group' in key:
                file_keys['group'].append(key)
            elif 'Subhalo' in key:
                file_keys['subhalo'].append(key)
            else:
                file_keys['part'].append(key)
        return file_keys

    def get_correct_keys(self, input_keys):
        corrected_keys = []
        for key in input_keys:
            if 'Group' in key:
                corrected_keys.append(f'Group/{key}')
            elif 'Subhalo' in key:
                corrected_keys.append(f'Subhalo/{key}')
            else:
                for i in self.part_types:
                    corrected_keys.append(f"PartType{i}/{key}")
        return corrected_keys
    
    #TODO
    def _check_input(self, path, part_types, keys):
        return
