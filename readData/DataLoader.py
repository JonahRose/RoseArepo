import readData.load_data as load
import numpy as np


class DataLoader():

    def __init__(self, path, part_types, keys):
        self._check_input(path, part_types, keys)

        self._input_keys = keys

        self.part_types = part_types
        self.path = path

        self.snap_path = ''
        self.group_path = ''
        self.get_paths()
        
        self.pt1_mass = self.get_pt1_mass()
        self.keys = [] #need to check what this does 

        #Change 'GroupMass' -> 'Groups/GroupMass'
        if type(part_types)==type([]):
            if len(part_types) > 0:
                self.keys = self.get_correct_keys()

        self.data = self._load_data()

    #TODO
    def __repr__(self):
        return ''

    #TODO: check if key is in data
    def __getattr__(self, attr):
        key = self._get_correct_keys([attr])
        return self.data[key]

    #TODO: have each just need .#.hdf5 so they can be passed to load.load_data()
    def get_paths(self):
        self.snap_path = ''
        self.group_path = ''
        return

    def get_pt1_mass(self):
        if 1 in part_types:
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

    def get_correct_keys(self):
        corrected_keys = []
        for key in self._input_keys:
            if 'Group' in key:
                corrected_keys.append(f'Group/{key}')
            elif 'Subhalo' in key:
                corrected_keys.append(f'Subhalo/{key}')
            else:
                for i in part_types:
                    corrected_keys.append(f"PartType{i}/{key}")
        return corrected_keys

    #TODO    
    def _load_data(self):
        return load.load_data(self.path, self.keys)
    
    #TODO
    def _check_input(self, path, part_types, keys):
        return
