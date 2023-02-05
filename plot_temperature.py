import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

example_trj = f'analysis/temperature/f_397_i_1200_mi_1.0_s_12_sgm_15_cut_21_c_10_t_100000_at_1000/bond_per_atoms_selected.out'
num_cycles = list(map(int, pd.read_csv(example_trj, header=None)[0].to_numpy() / 400000))

temperatures = [1000, 2000, 3000, 4000]


# Initialization
all_results = {}
for temp in temperatures:
    all_results[temp] = []  # each entry: [num, result]


# Load the dictionary
for temp in temperatures:
    for i, num in enumerate(num_cycles):
        entry = [num]
        trj_name = f'f_397_i_1200_mi_1.0_s_12_sgm_15_cut_21_c_10_t_100000_at_{temp}'
        data = pd.read_csv(f'analysis/temperature/{trj_name}/bond_per_atoms_selected.out', header=None)
        _, bonds_per_atom = data[0], data[1]
        entry.append(bonds_per_atom[i])
        all_results[temp].append(entry)
    all_results[temp] = np.array(all_results[temp])

plt.figure(figsize=(10, 8))
for prop, data in all_results.items():
    plt.plot(data[:, 0], data[:, 1], "-o", label=f'{prop} K', linewidth=2, markersize=6)
plt.xlabel(r"$N \rm _{cycle}$", fontsize=20)
plt.ylabel('bond per atom', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.show()
# all_results = utils.get_results(num_cycles, num_incs, seeds)
# all_results_mean_std = utils.get_mean_std(all_results)
# print(all_results)
# utils.plot_results(all_results_mean_std)