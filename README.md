# RoseArepo
This repo contains analysis and loading routines that can be used with outputs from Arepo simulations.

To import the data loading routines, you must import this class:

    from readData.DataLoader import DataLoader

Then the DataLoader can be called with these parameters:

    DataLoader(path, snap_num, part_types=-1, keys=[], sub_idx=-1, fof_idx=-1)

  path is the absolute or relative path to the folder containing the output data.

  part_types can be an integer or a list. If you do a list you have to specify the part type when you use particle data (cat['Coordinates'] vs cat['PartType4/Coordinates']).

  keys are the fields you want to load in as they are formatted here https://www.tng-project.org/data/docs/specifications/ .

  sub_idx is the subfind index.

  fof_idx is the group index.
  
  For part_types, sub_idx, and fof_idx, -1 means load everything 
 
If you want the positions and masses of all star particles in the first subhalo you would do:

    cat = DataLoader("path/to/output", 127, 4, ['Coordinates', 'Masses', 'SubhaloPos'], sub_idx=0)
    masses = cat['Masses']*1e10/cat.h
    pos = (cat['Coordinates'] - cat['SubhaloPos']) / cat.h
