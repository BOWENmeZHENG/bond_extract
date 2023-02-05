import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

example_trj = f'analysis/num_inc_high_sigma/f_397_i_400_mi_1.0_sgm_15_cut_21_c_10_t_100000_at_4000/bond_per_atoms_selected.out'
num_cycles = list(map(int, pd.read_csv(example_trj, header=None)[0].to_numpy() / 400000))

num_incs = [400, 600, 800, 1000, 1200]


# Initialization
all_results = {}
for num_inc in num_incs:
    all_results[num_inc] = []  # each entry: [num, result]


# Load the dictionary
for num_inc in num_incs:
    for i, num in enumerate(num_cycles):
        entry = [num]
        trj_name = f'f_397_i_{num_inc}_mi_1.0_sgm_15_cut_21_c_10_t_100000_at_4000'
        data = pd.read_csv(f'analysis/num_inc_high_sigma/{trj_name}/bond_per_atoms_selected.out', header=None)
        _, bonds_per_atom = data[0], data[1]
        entry.append(bonds_per_atom[i])
        all_results[num_inc].append(entry)
    all_results[num_inc] = np.array(all_results[num_inc])

plt.figure(figsize=(10, 8))
plt.text(0, 1.46, r'$\sigma=15.0 \rm \, \AA$', fontsize=20)
for prop, data in all_results.items():
    plt.plot(data[:, 0], data[:, 1], "-o", label=rf'$N_{{\rm inc}} = {prop}$', linewidth=2, markersize=6)
plt.xlabel(r"$N \rm _{cycle}$", fontsize=20)
plt.ylabel(r'$N \rm _{bond/atom}$', fontsize=20)
plt.ylim(top=1.47)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, frameon=False, loc='lower right')
plt.show()
# all_results = utils.get_results(num_cycles, num_incs, seeds)
# all_results_mean_std = utils.get_mean_std(all_results)
# print(all_results)
# utils.plot_results(all_results_mean_std)