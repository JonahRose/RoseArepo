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
            length -= boxsize
        if length < -boxlength / 2:
            length += boxsize

    return np.sqrt(dx*dx + dy*dy + dz*dz)

def get_box_cut(part_pos, center, length):
    xpr = center[0] + length
    xmr = center[0] - length
    ypr = center[1] + length
    ymr = center[1] - length
    zpr = center[2] + length
    zmr = center[2] - length

    is_x_in_box = (part_pos[:, 0] < xpr) & (part_pos[:, 0] > xmr)
    is_y_in_box = (part_pos[:, 1] < ypr) & (part_pos[:, 1] > ymr)
    is_z_in_box = (part_pos[:, 2] < zpr) & (part_pos[:, 2] > zmr)
    
    return is_x_in_box & is_y_in_box & is_z_in_box

def calc_dens(path, snap, gal_pos, min_r=0.1, max_r=800, r_step=1.05):

    cat = DataLoader(path, snap, part_types=[1,2], keys=['Coordinates','Masses'])
    
    part_pos = np.concatenate((cat['PartType2/Coordinates'], cat['PartType1/Coordinates']))
    part_mass = np.concatenate((cat['PartType2/Masses'], np.zeros(cat['PartType1/Coordinates'].shape[0])+cat.pt1_mass))

    #get the particles in the box 
    box_cut = get_box_cut(part_pos, gal_pos, max_r)
    part_pos  = part_pos[box_cut] - gal_pos
    part_mass = part_mass[box_cut] * 1e10/0.7
    part_r2 = np.sum(np.square(part_pos), axis=1)

    #calculate density
    dens_li = []
    all_r = []
    prev_r = 0
    r = min_r
    N = 4/3*np.pi
    while(r < max_r):
        r2 = r*r
        pr2 = prev_r*prev_r
        scut = (part_r2 < r2) & (part_r2 < pr2)
        in_s = part_pos[scut]
        in_s_mass = part_mass[scut]
        out_vol = N*r*r2
        in_vol = N*prev_r*pr2
        mass = np.sum(in_s_mass)
        
        dens_li.append(mass / (out_vol - in_vol))
        all_r.append(r)

        prev_r = r
        r *= r_step

    return all_r, dens_li

def calc_angular_momentum_vector(path, snap, gal_idx):
    #l = sum(r x p)
    
    gal_cat = DataLoader(path, snap, keys=["SubhaloPos", "SubhaloVel", "SubhaloHalfmassRad"])
    star_cat = DataLoader(path, snap, part_types=[4], keys=["Masses", "Coordinates", "Velocities"])

    pos = gal_cat['SubhaloPos'][gal_idx]
    vel = gal_cat["SubhaloVel"][gal_idx]
    max_rad = gal_cat["SubhaloHalfmassRad"][gal_idx]*2

    coords = star_cat["Coordinates"]
    pos_xcut = (coords[:,0] < pos[0]+max_rad) & (coords[:,0] > pos[0]-max_rad)
    pos_ycut = (coords[:,1] < pos[1]+max_rad) & (coords[:,1] > pos[1]-max_rad)
    pos_zcut = (coords[:,2] < pos[2]+max_rad) & (coords[:,2] > pos[2]-max_rad)
    pos_cut = pos_xcut & pos_ycut & pos_zcut

    star_pos = coords[pos_cut]
    star_vel = star_cat['Velocities'][pos_cut]
    star_mass = star_cat['Masses'][pos_cut]
    
    r_xyz = star_pos - pos
    v_xyz = star_vel - vel

    p_xyz = np.multiply(v_xyz, star_mass[:, np.newaxis])
    l_xyz = np.sum(np.cross(r_xyz, p_xyz), axis=0)

    return l_xyz

def get_next_gal(prev_mass, prev_loc, mass, pos, boxsize, tol=300):

    mcut = (mass < prev_mass*2) & (mass > prev_mass*0.1)

    idx_li = np.arange(mass.shape[0])[mcut] #just loop over galaxies with a reasonable mass
    if len(idx_li) == 0:
        return -1

    for idx in idx_li:  
        new_loc = pos[idx]
        dist = calc_dist(prev_loc, new_loc, boxsize)  
        if dist < tol:
            return idx
    if idx == idx_li[-1]:
        return -1


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
        if 'Subhalo/SubhaloMass' not in cat.data:
            return subhalo_idx_li
        pos = cat['SubhaloPos']
        mass = cat['SubhaloMass']
        mcut = (mass < prev_mass*2) & (mass > prev_mass*0.1)

        idx_li = np.arange(mass.shape[0])[mcut] #just loop over galaxies with a reasonable mass
        if len(idx_li) == 0:
            return subhalo_idx_li

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
