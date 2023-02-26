import numpy as np
import os
import utils

BL = 1.7

sigmas = [3, 5, 7, 9, 11, 13]
seeds = [13, 14, 15]

for seed in seeds:
    for sigma in sigmas:

        trj_name = f'S_3_f_198_i_200_mi_1.0_s_{seed}_sgm_{sigma}_cut_{sigma+3}_c_10_t_50000_at_2000'
        os.makedirs(f'analysis/S3/{trj_name}', exist_ok=True)

        steps, box_sizes, xyz = utils.parse(f'trajectories/distribution/S3/{trj_name}.lammpstrj')
        print(xyz.keys())
        np.savetxt(f'analysis/S3/{trj_name}/steps_all.out', steps, delimiter=',', fmt='%d')
        np.savetxt(f'analysis/S3/{trj_name}/box_sizes_all.out', np.array(list(box_sizes.values())), delimiter=',', fmt='%.5f')
        ids = np.arange(0, 2200000, 200000, dtype=int)


        bond_counts, num_atomss, bond_per_atoms = [], [], []
        for frame_id in ids:
            frame = xyz[frame_id]
            boxsize = box_sizes[frame_id]
            frame_name = f'analysis/S3/{trj_name}/xyz_{frame_id}.out'
            boxsize_name = f'analysis/S3/{trj_name}/box_size_{frame_id}.out'
            np.savetxt(boxsize_name, [boxsize])
            np.savetxt(frame_name, frame, delimiter=',')
            boxsize = np.loadtxt(boxsize_name)
            frame = np.loadtxt(frame_name, delimiter=',')

            bond_count, num_atoms, bond_per_atom = utils.comp_bond(frame, boxsize, BL)
            bond_counts.append(bond_count)
            num_atomss.append(num_atoms)
            bond_per_atoms.append(bond_per_atom)
            # print(bond_counts, num_atomss, bond_per_atoms)
        np.savetxt(f'analysis/S3/{trj_name}/bond_counts_selected.out', np.vstack((ids, bond_counts)).T, delimiter=',', fmt='%d')
        np.savetxt(f'analysis/S3/{trj_name}/bond_per_atoms_selected.out', np.vstack((ids, bond_per_atoms)).T, delimiter=',', fmt='%.5f')