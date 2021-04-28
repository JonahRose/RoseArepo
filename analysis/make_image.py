import matplotlib.pyplot as plt
import numpy as np
import RoseArepo.DataLoader as DataLoader
import analysis.analyze as analyze

def hist(path, snap_num, part_types, side_length=1, thickness=1, center=[0,0,0]):

    keys = ['Coordinates', 'Masses']
    cat = DataLoader(path, snap_num, part_types, keys)

    side_length *= cat.boxsize
    thickness *= cat.boxsize

    if type(center)==type([]):
        center = np.array(center)

    offset_pos = np.zeros((0,3))
    mass = np.zeros(0)
    for i in range(6):
        key = f'PartType{i}/Masses'
        if key in cat.data.keys():
            mass = np.concatenate((mass, cat[key]))
            offset_pos = np.concatenate((offset_pos, cat[f'PartType{i}/Coordinates']))

    is_in_box = analyze.get_box_cut(offset_pos, np.array([0,0,0]), side_length)
    in_box = offset_pos[is_in_box]
    mass_in_box = mass[is_in_box]

    is_in_thickness_1 = analyze.get_box_cut(offset_pos[is_in_box], np.array([0,0,0]), [side_length, side_length, thickness])
    is_in_thickness_2 = analyze.get_box_cut(offset_pos[is_in_box], np.array([0,0,0]), [side_length, thickness, side_length])
    is_in_thickness_3 = analyze.get_box_cut(offset_pos[is_in_box], np.array([0,0,0]), [thickness, side_length, side_length])

    xy_data = in_box[is_in_thickness_1]
    xz_data = in_box[is_in_thickness_2]
    yz_data = in_box[is_in_thickness_3]

    xy_weights = mass_in_box[is_in_thickness_1]
    xz_weights = mass_in_box[is_in_thickness_2]
    yz_weights = mass_in_box[is_in_thickness_3]

    return xy_data, xz_data, yz_data, xy_weights, xz_weights, yz_weights


