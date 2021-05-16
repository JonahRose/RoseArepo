import h5py
import numpy as np

class Offsets():

    def __init__(self, gpath, spath, sub_idx, fof_idx, num_gfiles, num_sfiles):
        
        if sub_idx == -1 and fof_idx == -1:
            return 

        self.grp_path = gpath
        self.snap_path = spath
        self.num_grp_files = num_gfiles
        self.num_snap_files = num_sfiles
        self.sub_idx = sub_idx
        self.fof_idx = fof_idx

        self.gal_file = None
        self.gal_before = None
        self.part_files = None
        self.particles_in_gal = None
        self.particles_before_gal = None
        self.get_gal_file_num()

        return 


    #get which file number our galaxy is on
    def get_gal_file_num(self):
        if self.sub_idx > -1:
            idx = self.sub_idx
            hkey = "Nsubgroups_ThisFile"
            dkey = "Subhalo/SubhaloLenType"
            akey = "Subhalo/SubhaloLen"
        elif self.fof_idx > -1:
            idx = self.fof_idx
            hkey = "Ngroups_ThisFile"
            dkey = "Group/GroupLenType"
            akey = "Group/GroupLen"
        else:
            return #not a single galaxy

        #get the number of each particle in the galaxy
        #get the number of particles in galaxies before it 
        #get which files each particle type covers 

        num_gals=0 #number of galaxies in files before idx
        gal_file=-1
        particles_before_gal = np.zeros(6)
        particles_in_gal = np.zeros(6)
        for i in range(self.num_grp_files):
            with h5py.File(f"{self.grp_path}{i}.hdf5", "r") as ofile:
                num_gals_this_file = ofile['Header'].attrs[hkey]
                if idx < num_gals+num_gals_this_file:
                    gal_file = i
                    particles_in_gal = np.array(ofile[dkey][idx-num_gals])
                    particles_before_gal_this_file = np.sum(ofile[dkey][:idx-num_gals], axis=0)
                    particles_before_gal += particles_before_gal_this_file
                    break
                if dkey not in ofile: #if a file has no data in it - found in dm only unifrom box
                    particles_this_file = 0
                else:
                    particles_this_file = np.array(ofile[dkey])
                particles_before_gal += np.sum(particles_this_file, axis=0)
                num_gals += num_gals_this_file

        self.gal_file = gal_file
        self.gal_before = num_gals

        self.particles_in_gal = particles_in_gal
        self.particles_before_gal = particles_before_gal

        particle_files = np.zeros(6)
        num_parts_before = np.zeros(6)
        num_parts_before_file = np.zeros(6)
        done = [False]*6
        for i in range(self.num_grp_files):
            if False not in done:
                break
            with h5py.File(f"{self.snap_path}{i}.hdf5", "r") as ofile:
                num_parts_file = np.array(ofile['Header'].attrs["NumPart_ThisFile"])
                for part in range(6):
                    if done[part]:
                        continue
                    if num_parts_file[part]+num_parts_before[part] >= particles_before_gal[part]+particles_in_gal[part]:
                        particle_files[part] = i
                        done[part] = True
                        num_parts_before_file[part] = particles_before_gal[part] - num_parts_before[part]
                    else:
                        num_parts_before[part] += num_parts_file[part]

        self.part_files = particle_files

        return

