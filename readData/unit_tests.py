from DataLoader import DataLoader
import numpy as np

def main():

    test_initialization()
    test_loading_non_output_data()
    test_offsets()

    return


def test_initialization():
    
    print("Testing Initialization.....", end='')
    cat = DataLoader('/home/j.rose/gr/CDM/', 127)
    print('\x1b[1;32;40m' + 'Passed' + '\x1b[0m')

    return

def test_loading_non_output_data():

    print("Testing Header Info.....", end='')
    cat = DataLoader('/home/j.rose/gr/CDM/', 127)

    h = cat.h == .6909
    boxsize = cat.boxsize == 100000
    time = cat.time == 0.9999999999999998 
    redshift = cat.redshift == 2.220446049250313e-16
    parts = np.sum(cat.num_parts == np.array([175744633,55717200,137129630,0,22991720,11])) == 6
    subs = cat.num_subhalos == 35066
    halos = cat.num_halos == 35066
    pfiles = cat.num_part_files == 64
    gfiles = cat.num_grp_files == 64

    pass_list = [h,boxsize,time,redshift,parts,subs,halos,pfiles,gfiles]

    if np.sum(pass_list) == 9:
        print('\x1b[1;32;40m' + 'Passed' + '\x1b[0m')
    else:
        print('\x1b[1;31;40m' + 'Failed' + '\x1b[0m')
    
    return

def test_offsets():

    print("Testing Particle Offsets.....", end='')

    diff = 0

    cat = DataLoader("/home/j.rose/gr/CDM/", 127, part_types=[1], fof_idx=0, sub_idx=0, keys=['ParticleIDs'])
    ids = cat['ParticleIDs'][:3].astype(int)
    correct = np.array([22252385, 25168921, 25860276])
    diff += np.sum(ids - correct)

    cat = DataLoader("/home/j.rose/gr/CDM/", 127, part_types=[1], fof_idx=10, sub_idx=10, keys=['ParticleIDs'])
    ids = cat['ParticleIDs'][:3].astype(int)
    correct = np.array([49527005, 49249142, 49387742])
    diff += np.sum(ids - correct)

    cat = DataLoader("/home/j.rose/gr/CDM/", 127, part_types=[4], fof_idx=30, sub_idx=0, keys=['ParticleIDs'])
    ids = cat['ParticleIDs'][:3].astype(int)
    correct = np.array([100283627291, 100383370441, 100272350342])
    diff += np.sum(ids - correct)

    cat = DataLoader("/home/j.rose/gr/CDM/", 127, part_types=[1], fof_idx=300, sub_idx=2, keys=['ParticleIDs'])
    ids = cat['ParticleIDs'][:3].astype(int)
    correct = np.array([47150913, 47150912, 47149262])
    diff += np.sum(ids - correct)

    if diff == 0:
        print('\x1b[1;32;40m' + 'Passed' + '\x1b[0m')
    else:
        print('\x1b[1;31;40m' + 'Failed' + '\x1b[0m')

    return


if __name__=="__main__":
    main()
