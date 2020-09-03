import matplotlib.pyplot as plt
import numpy as np
import RoseArepo.DataLoader as DataLoader

def hist(path, part_type):

    keys = []
    for pt in part_type:
        if pt != 1:
            keys.append(f'PartType{pt}/Masses')
        keys.append(f'PartType{pt}/Coordinates')
    Loader = DataLoader(path, part_type, keys)

    return
