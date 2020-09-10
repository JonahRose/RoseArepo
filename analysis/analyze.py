import matplotlib.pyplot as plt
import numpy as np
from readData.DataLoader import DataLoader

#accounts for wrapping around length of box
def calc_dist(loc1, loc2, boxlength):
    return

#returns the indecies for the largest halo in the box
#TODO: change to track any halo
#TODO: finish
def track_largest(path, start_snap, end_snap=0, tol=300):

    keys = ['SubhaloPos', 'SubhaloMass']
    cat_0 = DataLoader(path, start_snap, keys=keys)

    subhalo_idx_li = [0]
    prev_mass = cat_0['SubhaloMass'][0]
    prev_loc = cat_0['SubhaloPos'][0]

    for snap_num in range(start_snap-1, end_snap-1, -1):
        cat = DataLoader(path, snap_num, keys)
        new_loc = cat['SubhaloPos'][0]
        dist = calc_dist(prev_loc, new_loc, cat.boxsize) #TODO: get boxsize 

        if dist < tol:
            subhalo_idx_li.append(idx)

    return subhalo_idx_li
