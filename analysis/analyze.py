import matplotlib.pyplot as plt
import numpy as np
from readData.DataLoader import DataLoader
import ctypes

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

    if type(length) == type([]):
        xlen = length[0]
        ylen = length[1]
        zlen = length[2]
    else:
        xlen = ylen = zlen = length

    xpr = center[0] + xlen
    xmr = center[0] - xlen
    ypr = center[1] + ylen
    ymr = center[1] - ylen
    zpr = center[2] + zlen
    zmr = center[2] - zlen

    is_x_in_box = (part_pos[:, 0] < xpr) & (part_pos[:, 0] > xmr)
    is_y_in_box = (part_pos[:, 1] < ypr) & (part_pos[:, 1] > ymr)
    is_z_in_box = (part_pos[:, 2] < zpr) & (part_pos[:, 2] > zmr)
    
    return is_x_in_box & is_y_in_box & is_z_in_box

def checklen(x):
    return len(np.array(x,ndmin=1));
def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (abs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (abs(input)<=xmax);

def int_round(x):
    return np.int(np.round(x));

def get_particle_hsml(x, y, z, DesNgb=32, Hmax=0):
    x = fcor(x)
    y = fcor(y)
    z = fcor(z)
    N = checklen(x)

    ok = ok_scan(x) & ok_scan(y) & ok_scan(z)
    x = x[ok]
    y = y[ok]
    z = z[ok]

    if Hmax == 0:
        dx = np.max(x) - np.min(x)
        dy = np.max(y) - np.min(y)
        dz = np.max(z) - np.min(z)
        ddx = np.max([dx, dy, dz])
        Hmax = 5*ddx*np.power(np.float(N), -1/3)

    exec_call = '/home/j.rose/torreylabtools_git/C/c_libraries/' + '/StellarHsml/starhsml.so'
    h_routine = ctypes.cdll[exec_call]

    h_out_cast = ctypes.c_float*N
    H_OUT = h_out_cast()

    h_routine.stellarhsml(ctypes.c_int(N),
                          vfloat(x), vfloat(y), vfloat(z),
                          ctypes.c_int(DesNgb),
                          ctypes.c_float(Hmax),
                          ctypes.byref(H_OUT) )
        
    h = np.ctypeslib.as_array(H_OUT)
    return h


def calc_dens(path, snap, gal_pos, min_r=0.1, max_r=800, r_step=1.05):

    cat = DataLoader(path, snap, part_types=[1,2], keys=['Coordinates','Masses'])
    
    part_pos = np.concatenate((cat['PartType2/Coordinates'], cat['PartType1/Coordinates']))
    part_mass = np.concatenate((cat['PartType2/Masses'], np.zeros(cat['PartType1/Coordinates'].shape[0])+cat.pt1_mass))

    #get the particles in the box 
    box_cut = get_box_cut(part_pos, gal_pos, max_r)
    part_pos  = part_pos[box_cut] - gal_pos
    part_mass = part_mass[box_cut] * 1e10/cat.h
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
        scut = (part_r2 < r2) & (part_r2 > pr2)
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

def calc_surface_dens_profile(path, snap, gal_idx, min_r, max_r, r_step=1.05, part_types=[4], zheight=None):
    cat = DataLoader(path, snap, part_types=part_types, keys=['Coordinates','Masses', 'Velocities'])
    gal_cat = DataLoader(path, snap, keys=['SubhaloPos', 'SubhaloVel'])

    gal_pos = gal_cat['SubhaloPos'][gal_idx]/cat.h
    gal_vel = gal_cat['SubhaloVel'][gal_idx]

    pos = cat['Coordinates']/cat.h
    mass = cat['Masses']*1e10/cat.h
    vel = cat['Velocities']

    box_cut = get_box_cut(pos, gal_pos, max_r)   
    pos_box = pos[box_cut] - gal_pos
    vel_box = vel[box_cut] - gal_vel
    mass_box = mass[box_cut]

    part_r2 = np.sum(np.square(pos_box), axis=1)

    pos_box, vel_box = get_rotate_data(pos_box, vel_box, mass_box, face_on=True)

    #l_hat = calc_angular_momentum_vector(path, snap, 0)
    #pos_box = rotate_data(pos_box, l_hat)
   

    if zheight is not None:
        zcut = (pos_box[:,2] > -zheight) & (pos_box[:,2] < zheight)
        pos_box = pos_box[zcut]
        mass_box = mass_box[zcut]
        part_r2 = np.sum(np.square(pos_box), axis=1)

    dens_li = []
    all_r = []
    prev_r = 0
    r = min_r
    while(r < max_r):
        r2 = r*r
        pr2 = prev_r*prev_r
        scut = (part_r2<r2) & (part_r2 > pr2)
        in_s = pos_box[scut]
        in_s_mass = mass_box[scut]
        out_area = np.pi*r2
        in_area = np.pi*pr2
        mass = np.sum(in_s_mass)

        dens_li.append(mass / (out_area - in_area))
        all_r.append(r)

        prev_r = r
        r *= r_step

    return all_r, dens_li

#adapted from torreylabtools
#assumes data has already been centered and had galaxy's velocity removed
def get_rotate_data(coords, velocities, masses, phi=0, theta=0, edge_on=False, face_on=False, get_pt=False, r_max=5, r_min=0):

    x = coords[:,0]
    y = coords[:,1]
    z = coords[:,2]
    vx = velocities[:,0]
    vy = velocities[:,1]
    vz = velocities[:,2]

    if edge_on or face_on:
        
        r2 = x*x+y*y+z*z
        scut = (r2<r_max*r_max) & (r2>r_min*r_min)

        lz = np.sum(masses[scut] * (x[scut]*vy[scut] - y[scut]*vx[scut]))
        lx = np.sum(masses[scut] * (y[scut]*vz[scut] - z[scut]*vy[scut]))
        ly = np.sum(masses[scut] * (z[scut]*vx[scut] - x[scut]*vz[scut]))

        phi = np.arctan2(ly, lx)
        theta = np.arctan2(np.sqrt(lx*lx + ly*ly), lz)

        if edge_on:
            #phi += np.pi/2
            theta += np.pi/2

    if get_pt:
        return phi, theta

    x_ = -z  * np.sin(theta) + (x * np.cos(phi) + y *np.sin(phi)) * np.cos(theta)
    y_ = -x  * np.sin(phi)   + y  * np.cos(phi)
    z_ =  z  * np.cos(theta) + (x * np.cos(phi) + y *np.sin(phi)) * np.sin(theta)
    vx_ = -vz  * np.sin(theta) + (vx * np.cos(phi) + vy *np.sin(phi)) * np.cos(theta)
    vy_ = -vx  * np.sin(phi)   + vy  * np.cos(phi)
    vz_ =  vz  * np.cos(theta) + (vx * np.cos(phi) + vy *np.sin(phi)) * np.sin(theta)

    coords[:,0] = x_
    coords[:,1] = y_
    coords[:,2] = z_

    velocities[:,0] = vx_
    velocities[:,1] = vy_
    velocities[:,2] = vz_

    return coords, velocities

def get_next_gal(prev_mass, prev_loc, mass, pos, boxsize, tol=300):

    mcut = (mass < prev_mass*2) & (mass > prev_mass*0.5)

    idx_li = np.arange(mass.shape[0])[mcut] #just loop over galaxies with a reasonable mass
    if len(idx_li) == 0:
        return -1

    all_dist = []
    for idx in idx_li:  
        new_loc = pos[idx]
        dist = calc_dist(prev_loc, new_loc, boxsize)  
        all_dist.append(dist)
        #if dist < tol:
        #    return idx
    all_dist = np.array(all_dist)
    if min(all_dist) < tol:
        return idx_li[all_dist==min(all_dist)][0]
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
