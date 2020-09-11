from readData.DataLoader import DataLoader
import analysis.analyze as analyze

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def main():

    path = '/home/j.rose/Projects/my_MW_DAO/IllustrisTNGcode/M2812/'
    start_snap = 126
    end_snap = 0
    gal_idx = 0

    plot_sfr = True
    plot_mass_type = True
    plot_dens = [126, 110, 90, 70, 50]
    plot_pos = True

    target_idx = analyze.track_largest(path, start_snap, end_snap, start_gal_idx=gal_idx, tol=600)
    x_data = {'dens':[]}
    y_data = {'sfr':[], 'mass_type':[], 'dens':[], 'pos':[]}
    time = []
    dens_time = []
    idx_counter = 0

    for i in range(0, len(target_idx)):
        snap = start_snap - i
        if snap % 10 == 0:
            print(snap)
        idx = target_idx[idx_counter]
        if idx == -1:
            continue
        idx_counter += 1

        keys = []
        if plot_sfr:
            keys.append('SubhaloSFR')
        if plot_mass_type:
            keys.append('SubhaloMassType')
        if plot_dens is not False or plot_pos:
            keys.append('SubhaloPos')

        cat = DataLoader(path, snap, keys=keys)
        time.append(cat.time)

        if plot_sfr:
            y_data['sfr'].append(cat['SubhaloSFR'][idx])
        if plot_mass_type:
            y_data['mass_type'].append(cat['SubhaloMassType'][idx])
        if plot_dens is not False:
            if snap in plot_dens:
                pos = cat['SubhaloPos'][idx]
                dens_time.append(cat.time)
                r, d = analyze.calc_dens(path, snap, pos)
                x_data['dens'].append(r)
                y_data['dens'].append(d)
        if plot_pos:
            y_data['pos'].append(cat['SubhaloPos'][idx])

    for key in y_data.keys():
        y_data[key] = np.array(y_data[key])
    for key in x_data.keys():
        x_data[key] = np.array(x_data[key])

    num_frames = np.sum([plot_sfr, plot_mass_type, plot_pos])
    num_frames += 0 if plot_dens is False else 1
    fig, ax = plt.subplots(num_frames, figsize=(10, 5*num_frames))
    ax_idx = 0

    if plot_sfr:
        ax[ax_idx].plot(time, y_data['sfr'], 'b-')
        ax[ax_idx].set_title("SFR")
        ax[ax_idx].set_xscale("log")
        ax[ax_idx].set_xlabel("Time")
        ax[ax_idx].set_ylabel("SFR")
        ax_idx += 1

    if plot_mass_type:
        for i in range(6):
            ax[ax_idx].plot(time, y_data['mass_type'][:,i], label=f'PartType{i}')
        ax[ax_idx].legend()
        ax[ax_idx].set_title("Mass Breakdown")
        ax[ax_idx].set_xscale("log")
        ax[ax_idx].set_yscale("log")
        ax[ax_idx].set_xlabel("Time")
        ax[ax_idx].set_ylabel("Mass")
        ax_idx += 1       

    if plot_dens is not False:
        for i in range(len(x_data['dens'])):
            ax[ax_idx].plot(x_data['dens'][i], y_data['dens'][i], label=f"a = {round(dens_time[i],3)}")
        ax[ax_idx].legend()
        ax[ax_idx].set_title("Radial Density")
        ax[ax_idx].set_xscale("log")
        ax[ax_idx].set_yscale("log")
        ax[ax_idx].set_xlabel("Radius (kpc)")
        ax[ax_idx].set_ylabel("Density ($M_\odot$ / $kpc^3$)")
        ax_idx += 1       

    if plot_pos:
        colors = plt.cm.get_cmap('inferno')
        scat = ax[ax_idx].scatter(y_data['pos'][:,0], y_data['pos'][:, 1], c=time, cmap=colors, norm=LogNorm())
        plt.colorbar(scat)
        ax[ax_idx].set_xlabel("x Position")
        ax[ax_idx].set_ylabel("y Position")
        ax[ax_idx].set_title("Position")
        ax_idx += 1
    fig.savefig("plots/tracker.pdf")

    return

if __name__=="__main__":
    main()
