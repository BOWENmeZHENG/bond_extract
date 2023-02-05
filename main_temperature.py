import numpy as np
import os
import utils

BL = 1.7

temperatures = [1000, 2000, 3000, 4000]

for temp in temperatures:

        trj_name = f'f_397_i_1200_mi_1.0_s_12_sgm_15_cut_21_c_10_t_100000_at_{temp}'
        os.makedirs(f'analysis/temperature/{trj_name}', exist_ok=True)

        steps, box_sizes, xyz = utils.parse(f'trajectories/temperature/{trj_name}.lammpstrj', line_limit=13000000)
        np.savetxt(f'analysis/temperature/{trj_name}/steps_all.out', steps, delimiter=',', fmt='%d')
        np.savetxt(f'analysis/temperature/{trj_name}/box_sizes_all.out', np.array(list(box_sizes.values())), delimiter=',', fmt='%.5f')
        ids = np.arange(0, 4400000, 400000, dtype=int)


        bond_counts, num_atomss, bond_per_atoms = [], [], []
        for frame_id in ids:
            frame = xyz[frame_id]
            boxsize = box_sizes[frame_id]
            frame_name = f'analysis/temperature/{trj_name}/xyz_{frame_id}.out'
            boxsize_name = f'analysis/temperature/{trj_name}/box_size_{frame_id}.out'
            np.savetxt(boxsize_name, [boxsize])
            np.savetxt(frame_name, frame, delimiter=',')
            boxsize = np.loadtxt(boxsize_name)
            frame = np.loadtxt(frame_name, delimiter=',')

            bond_count, num_atoms, bond_per_atom = utils.comp_bond(frame, boxsize, BL)
            bond_counts.append(bond_count)
            num_atomss.append(num_atoms)
            bond_per_atoms.append(bond_per_atom)
            print(bond_counts, num_atomss, bond_per_atoms)
        np.savetxt(f'analysis/temperature/{trj_name}/bond_counts_selected.out', np.vstack((ids, bond_counts)).T, delimiter=',', fmt='%d')
        np.savetxt(f'analysis/temperature/{trj_name}/bond_per_atoms_selected.out', np.vstack((ids, bond_per_atoms)).T, delimiter=',', fmt='%.5f')