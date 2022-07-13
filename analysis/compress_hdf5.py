from mpi4py import MPI
import numpy as np
import os,h5py,time
import functools

"""
This script will copy an existing snapshot file and write a new one in a compressed format.
The script is modified from a similar script provided by Francisco Antonio Villaescusa Navarro.

With the default compression value, users should expect a ~60% reduction in the size of the file
    and a ~3x increase in readtime.
Users can use 'h5diff' in the command line (or set check_files to True) to ensure that the 
    compressed files are not different from the originals

This script can take a while to execute, it is recomended to run in parallel with processors
    roughly equal to the number of individual files in a single snapshot

The script can be run with something similar to from the command line:
    mpirun -n 32 python compress_hdf5.py
"""
def main():

    ###### MPI DEFINITIONS ######                                    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    ###################################### INPUT ###########################################
    root        = './' #path to snapshots
    gzip        = 4 #4 is default, 0 is no compression, 9 is max compression
    snapnums    = np.arange(0,128) #which snapshots to compress (will skip if not present)
    check_files = False #run h5diff after the compressed file is written

    folder_out  = 'compressed/' #compressed snapshot will be written to {root}/{folder_out}/
    num_files   = 32 #how many files were written in parallel when the snapshot was written
    folder_base = 'snapdir' #base name for the snapshot directory
    file_base   = 'snap' #base name for the individual snapshot files
    ########################################################################################

    #make the output folder if it does not already exist
    if rank == 0:
        if not os.path.exists(f"{root}/{folder_out}"):
            os.mkdir(f"{root}/{folder_out}")
    else:
        time.sleep(1)

    #spread the jobs over the different tasks
    realizations, snapnums = get_realizations(num_files, snapnums, size, rank)
    
    # do a loop over the different snapshots
    for snap in snapnums:

        #get paths to current snapshot
        input_base, output_base, file_name = get_paths(root, folder_out, num_files, 
                                                        folder_base, file_base, snap)

        #flag = 0: do nothing, flag = 1: continue to next snapshot
        flag = do_folder_checks(input_base, output_base, num_files, snap, rank)
        if flag == 1:
            continue

        # do a loop over inividual files in snapshot folder (if it is a folder)
        for file_num in realizations:

            start = time.time()

            # find the names of the original and compressed snapshot
            snap_in  = input_base + file_name
            snap_out = output_base + file_name

            if num_files > 1:
                snap_in  += f'{file_num}.hdf5'
                snap_out += f'{file_num}.hdf5'

            if do_file_checks(root, folder_out, snap_in, snap_out):
                continue

            # compress the snapshot
            compress_snapshot(snap_in, snap_out, gzip, check_files)

            print(f'{rank}: Compressing {snap_in} with gzip={gzip} \n' \
                  f'    Time taken = {time.time() - start:.3f}')

    return

#get the path names for input and output data
def get_paths(root, folder_out, num_files, folder_base, file_base, snap):

    #get the paths to input and output files
    if num_files == 1:
        input_base  = f"{root}/"
        output_base = f"{root}/{folder_out}/"
        file_name   = f"{file_base}_{snap:03d}.hdf5"
    elif num_files > 1:
        input_base  = f"{root}/{folder_base}_{snap:03d}/"
        output_base = f"{root}/{folder_out}/{folder_base}_{snap:03d}/"
        file_name   = f"{file_base}_{snap:03d}." 
    else:
        raise ValueError(f"Number of Files: {num_files} not recognized")

    return input_base, output_base, file_name

#check that the input and output folders are present
def do_folder_checks(input_base, output_base, num_files, snap, rank):

    #if the input snapshot is not present, skip it
    if not os.path.exists(input_base):
        #cannot find folder with single-file snapshots
        if str(snap) not in input_base:
            raise Exception(f"Cannot find {input_base}")
        #input snapshot not present, continue to next
        if rank == 0:
            print(f"Did not find Snapshot {snap}")
        return 1

    #if the output snapshot folder does not exist, create it
    if num_files > 1:
        if not os.path.exists(output_base):
            os.mkdir(output_base)

    return 0

#get the current file for this task
def get_realizations(num_files, snapnums, size, rank):

    if num_files == 1:
        numbers      = snapnums
        indexes      = np.where(np.arange(numbers.shape[0])%size==rank)[0]
        snapnums     = numbers[indexes]
        realizations = [1]
    elif num_files > 1:
        numbers      = np.arange(num_files)
        indexes      = np.where(np.arange(numbers.shape[0])%size==rank)[0]
        realizations = numbers[indexes]
    else:
        raise Exception("Unknown number of files")

    return realizations, snapnums

#checks if the file is already present and can be written
def do_file_checks(root, folder_out, snap_in, snap_out):

    if not os.path.exists(snap_in):
        return True

    if os.path.exists(snap_out):  
        print(f"File {snap_out} already present")
        return True

    if not os.path.exists(root+folder_out):
        raise Exception(f"Output folder {root+folder_out} does not exist")

    return False

# This function takes a snapshot and compress it
def compress_snapshot(snap_in, snap_out, gzip, check_files): 

    # open the input and output files
    f_in  = h5py.File(snap_in,  'r')
    f_out = h5py.File(snap_out, 'w')

    # core call
    copy_datasets = functools.partial(copy_level0, f_in, f_out, gzip)
    f_in.visititems(copy_datasets)
    f_in.close();  f_out.close()

    # check if input and output files are the same
    if check_files:  
        os.system(f'h5diff {snap_in} {snap_out}')

    return

# This function writes the compressed file
def copy_level0(fs, fd, compression_opts, name, node):
    """full copy of datasets + header: chunking + compressing"""

    if isinstance(node, h5py.Dataset):
        compression = 'gzip'
        if not compression_opts:
            compression = None

        dnew = fd.create_dataset(
            name, data=node, dtype=node.dtype, chunks=True,
            shuffle=True, compression=compression,
            compression_opts=compression_opts, fletcher32=True)

    elif isinstance(node, h5py.Group) and name=='Header':
        fs.copy(name, fd, name=name)

    #I'm not sure what falls under else, but it doesn't seem to cause problems
    else:
        pass 
        #print("Group", name)

    return

if __name__=="__main__":
    main()
