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
        self.sub_idx = int(sub_idx)
        self.fof_idx = int(fof_idx)

        self.gal_file = None
        self.gal_before = None
        self.part_files = None
        self.particles_in_gal = None
        self.particles_before_gal = None
        self.get_gal_file_num()

        return 


    #get which file number our galaxy is on and the particle locations
    #probably not work if there are subhalos that are not in groups
    #may not work if there are fofs without subhalos before fofs with subhalos 
    def get_gal_file_num(self):
        dkey = "Group/GroupLenType"
        if self.sub_idx > -1:
            idx = self.sub_idx
            pkey = "Subhalo/SubhaloLenType"
            gkey = "Subhalo/SubhaloGrNr"
            is_subs = True
            hkey = "Nsubgroups_ThisFile"
        elif self.fof_idx > -1:
            idx = self.fof_idx
            pkey = "Group/GroupLenType"
            is_subs = False
            hkey = "Ngroups_ThisFile"
        else:
            return #not a single galaxy

        #get the number of each particle in the galaxy
        #get the number of particles in galaxies before it 
        #get which files each particle type covers 

        num_groups = 0
        num_gals=0 #number of subs/fofs in files before idx
        num_gals_groups = 0 #number of subhalos before our fof, using the fof groups to count
        gal_file=-1 #the file that the sub/fof is located on
        particles_before_gal = np.zeros(6) #number of particles of each type before the galaxy
        particles_in_gal = np.zeros(6) #number of particles on each type in the galaxy
        subs_still_to_go = None
        for i in range(self.num_grp_files):
            with h5py.File(f"{self.grp_path}{i}.hdf5", "r") as ofile:
                #get the number of subs/fofs in this file
                num_gals_this_file = int(ofile['Header'].attrs[hkey])
                num_groups_this_file = int(ofile['Header'].attrs['Ngroups_ThisFile'])
                #print(num_groups_this_file)
                if num_gals_this_file == 0: #if a file has no data in it - found in dm only unifrom box
                    continue
               
                #the number of subhalos that are already accounted for in this file
                #this will happen when the fof that a subhalo is in is not the first on the file 
                num_subs_this_file = 0

                #if the groups for the subs on this file have already been accounted for, then remove
                #the galaxies from the counters that are relative to the current file index 
                if subs_still_to_go is not None:
                    #print("test", subs_still_to_go, num_gals_this_file)
                    if subs_still_to_go >= num_gals_this_file:
                        subs_before_group -= num_gals_this_file
                        subs_still_to_go -= num_gals_this_file
                        num_gals += num_gals_this_file
                        continue

                #print("test1", idx, num_gals, num_gals_this_file, subs_still_to_go, particles_before_gal)
                #print(ofile['Group/GroupFirstSub'][-1], ofile['Subhalo/SubhaloGrNr'][-1])
                if idx < num_gals+num_gals_this_file: #if the sub/fof is on this file
                    gal_file = i
                    #print(i, idx, num_gals, num_gals_this_file)
                    
                    #get any particles on this file in fof groups before out sub's group
                    if is_subs:

                        #if we have had to get the number of subhalos before ours from a different file
                        if subs_still_to_go is not None:
                            #print("test2", subs_before_group, subs_still_to_go)
                            particles_before_gal += np.sum(ofile['Subhalo/SubhaloLenType'][subs_before_group:subs_still_to_go], axis=0)
                            particles_in_gal = np.array(ofile[pkey][idx-num_gals])
                            break
                        #if the subhalo is on the same file as its fof
                        #account for any groups that are before our subhalo
                        else:
                            #print("test")
                            gidx = int(ofile[gkey][0]) #num_groups# int(ofile[gkey][0])
                            group_counter = 0
                            #print(idx, num_gals, gidx, ofile[gkey][idx-num_gals])
                            #print(i, num_groups, list(ofile['Group/GroupFirstSub']),list(ofile['Group/GroupNsubs']))
                            #if gidx > 0:
                            #    num_subs_this_file += ofile['Group/GroupFirstSub'][gidx-1]
                            #    num_subs_this_file += ofile['Group/GroupNsubs'][gidx-1]
                            #print(gidx, idx, num_gals, num_gals_this_file, num_groups_this_file, ofile[gkey][idx-num_gals], ofile[gkey][0], num_groups, num_gals_groups)
                            #print(ofile['Group/GroupFirstSub'][0], ofile['Subhalo/SubhaloGrNr'][0], ofile['Subhalo/SubhaloGrNr'][idx-num_gals], gidx, gidx-num_groups)
                            num_subs_this_file = ofile['Group/GroupFirstSub'][ofile[gkey][idx-num_gals]-num_groups]-num_gals
                            #print(subs_this_file, num_groups, ofile['Subhalo/SubhaloGrNr'][idx-num_gals])
                            
                            #while gidx < ofile[gkey][idx-num_gals]:
                            while num_groups < ofile['Subhalo/SubhaloGrNr'][idx-num_gals]:
                                #print("test")
                                particles_before_gal += np.array(ofile['Group/GroupLenType'][group_counter])
                                #print(particles_before_gal)
                                num_subs_this_group = int(ofile['Group/GroupNsubs'][group_counter])
                                #if num_subs_this_group != -1:
                                #    num_subs_this_file += num_subs_this_group
                                num_groups += 1
                                group_counter += 1


                    #get the number of particles in the current fof before our sub and the particles in
                    #the current fof/sub
                    particles_in_gal = np.array(ofile[pkey][idx-num_gals])
                    #print(particles_before_gal)
                    particles_before_gal += np.sum(ofile[pkey][num_subs_this_file:idx-num_gals], axis=0)
                    #print("test6", particles_before_gal, num_subs_this_file, idx-num_gals,  np.sum(ofile[pkey][num_subs_this_file:idx-num_gals], axis=0))
                    break

                #check if the group starts on this file but the subhalo is on a different one
                if is_subs:
                    nsubs = np.array(ofile['Group/GroupNsubs'])
                    #print(idx, num_gals, np.sum(nsubs), num_gals_groups)
                    if idx < num_gals+np.sum(nsubs) or (num_gals_groups >= num_gals and idx < num_gals_groups + np.sum(nsubs)):
                        #get the index into the subhalo's fof group
                        first_sub = np.array(ofile['Group/GroupFirstSub'])
                        group_li = np.arange(ofile['Header'].attrs['Ngroups_ThisFile'])
                        #print(group_li, num_gals_groups, idx, first_sub)
                        my_group_cut = (first_sub <= idx) & (first_sub != -1)
                        my_group = group_li[my_group_cut][-1]
                        #print("group", my_group, num_groups_this_file, first_sub[my_group], idx)

                        #find the number of particles and subhalos before the target one
                        particles_before_gal += np.sum(ofile['Group/GroupLenType'][:my_group], axis=0)
                        #print("test3", my_group, first_sub[my_group], nsubs[my_group])
                        #print(np.sum(nsubs[:my_group]) , num_gals_groups ,num_gals , num_gals_this_file) 
                        subs_before_group = np.sum(nsubs[:my_group]) + num_gals_groups -num_gals - num_gals_this_file 
                        #if subs_before_group < 0:
                        #    subs_before_group += num_gals_this_file - np.sum(nsubs[:my_group])
        
                        #figure out how many subs are on this file and subsequent files that we 
                        #still need to account for
                        subs_still_to_go = idx - num_gals - num_gals_this_file
                        subs_in_file = num_gals_this_file - (first_sub[my_group] - num_gals)
                        #print("test4", my_group, subs_before_group, subs_still_to_go, subs_in_file)

                        #subhalos before ours that are in the same fof but on a different file
                        if subs_in_file <=0: 
                            pass 
                        else:
                            last_parts = np.sum(ofile['Subhalo/SubhaloLenType'][-subs_in_file:], axis=0)
                            particles_before_gal += last_parts

                        num_gals += num_gals_this_file
                        continue
                
                #if the fof is not on this file add all particles on the file (parts in all fofs)
                particles_this_file = np.array(ofile['Group/GroupLenType']) #particles in groups
                #print("test5", particles_before_gal, np.sum(particles_this_file, axis=0))
                #print("len", len(particles_this_file), particles_before_gal)
                particles_before_gal += np.sum(particles_this_file, axis=0)
                #print("test5", particles_before_gal)
                num_gals += num_gals_this_file
                num_gals_groups += np.sum(ofile['Group/GroupNsubs']) 
                num_groups += num_groups_this_file

        #print(particles_in_gal, particles_before_gal)

        self.gal_file = gal_file
        self.gal_before = num_gals

        #print(particles_in_gal, particles_before_gal)
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

