import matplotlib.pyplot as plt
import numpy as np
from readData.DataLoader import DataLoader

#accounts for wrapping around length of box
def calc_dist(loc1, loc2, boxlength):
    dx = loc1[0] - loc2[0]
    dy = loc1[1] - loc2[1]
    dz = loc1[2] - loc2[2]

    for length in [dx, dy, dz]:
        if length > boxlength / 2:
            length -= boxdize
        if length < -boxlength / 2:
            length += boxdize

    return np.sqrt(dx*dx + dy*dy + dz*dz)

def calc_dens():
    return

#returns the indecies for the largest subhalo in the box
#returns -1 if a halo is not found within tolerance 
#TODO: change to handle mergers better
#          ie. make it so that if there are multiple in tol, choose best
def track_largest(path, start_snap, end_snap=0, start_gal_idx=0, tol=300):

    keys = ['SubhaloPos', 'SubhaloMass']
    cat_0 = DataLoader(path, start_snap, keys=keys)

    subhalo_idx_li = [start_gal_idx]
    prev_mass = cat_0['SubhaloMass'][start_gal_idx]
    prev_loc = cat_0['SubhaloPos'][start_gal_idx]
    for snap_num in range(start_snap-1, end_snap-1, -1):
        cat = DataLoader(path, snap_num, keys=keys)
        pos = cat['SubhaloPos']
        mass = cat['SubhaloMass']
        mcut = (mass < prev_mass*2) & (mass > prev_mass*0.5)

        idx_li = np.arange(mass.shape[0])[mcut] #just loop over galaxies with a reasonable mass
        for idx in idx_li:  
            new_loc = pos[idx]
            dist = calc_dist(prev_loc, new_loc, cat.boxsize)  
            if dist < tol:
                prev_loc = new_loc
                subhalo_idx_li.append(idx)
                break
        if idx == idx_li[-1]:
            subhalo_idx_li.append(-1)
                
    return subhalo_idx_li
