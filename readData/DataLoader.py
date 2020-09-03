import readData.load_data as load
import numpy as np


class DataLoader():

    def __init__(self, path, part_types, keys):
        self._check_input(path, part_types, keys)

        self.path = path
        self.input_keys = keys
        self.pt1_mass = 0

        if 1 in part_types:
            load.get_pt1_mass(self.path)

    def get_correct_keys(self):
        for key in self.input_keys:
            pass
        return

    def _load_data(self):
        return

    def _check_input(self, path, part_types, keys):
        return
