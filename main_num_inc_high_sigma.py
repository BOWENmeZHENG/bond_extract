import numpy as np
import os
import utils

BL = 1.7

num_incs = [400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
task = 'num_inc_high_sigma'

for num in num_incs:

        trj_name = f'f_397_i_{num}_mi_1.0_sgm_15_cut_21_c_10_t_100000_at_4000'
        os.makedirs(f'analysis/{task}/{trj_name}', exist_ok=True)

        steps, box_sizes, xyz = utils.parse(f'trajectories/{task}/{trj_name}.lammpstrj', line_limit=23000000)
        np.savetxt(f'analysis/{task}/{trj_name}/steps_all.out', steps, delimiter=',', fmt='%d')
        np.savetxt(f'analysis/{task}/{trj_name}/box_sizes_all.out', np.array(list(box_sizes.values())), delimiter=',', fmt='%.5f')
        ids = np.arange(0, 4400000, 400000, dtype=int)


        bond_counts, num_atomss, bond_per_atoms = [], [], []
        for frame_id in ids:
            frame = xyz[frame_id]
            boxsize = box_sizes[frame_id]
            frame_name = f'analysis/{task}/{trj_name}/xyz_{frame_id}.out'
            boxsize_name = f'analysis/{task}/{trj_name}/box_size_{frame_id}.out'
            np.savetxt(boxsize_name, [boxsize])
            np.savetxt(frame_name, frame, delimiter=',')
            boxsize = np.loadtxt(boxsize_name)
            frame = np.loadtxt(frame_name, delimiter=',')

            bond_count, num_atoms, bond_per_atom = utils.comp_bond(frame, boxsize, BL)
            bond_counts.append(bond_count)
            num_atomss.append(num_atoms)
            bond_per_atoms.append(bond_per_atom)
            print(bond_counts, num_atomss, bond_per_atoms)
        np.savetxt(f'analysis/{task}/{trj_name}/bond_counts_selected.out', np.vstack((ids, bond_counts)).T, delimiter=',', fmt='%d')
        np.savetxt(f'analysis/{task}/{trj_name}/bond_per_atoms_selected.out', np.vstack((ids, bond_per_atoms)).T, delimiter=',', fmt='%.5f')