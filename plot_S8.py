import utils
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_results(num_cycles, property_list, seed_list):
    # Initialization
    all_results = {}
    for prop in property_list:
        all_results[prop] = []  # each entry: [num, [results for all seeds]]


    # Load the dictionary
    for prop in property_list:
        for i, num in enumerate(num_cycles):
            entry = [num, []]
            for seed in seed_list:
                trj_name = f'S_8_f_195_i_200_mi_1.0_s_{seed}_sgm_{prop}_cut_{prop+3}_c_10_t_50000_at_2000'
                data = pd.read_csv(f'analysis/S8/{trj_name}/bond_per_atoms_selected.out', header=None)
                _, bonds_per_atom = data[0], data[1]
                entry[-1].append(bonds_per_atom[i])
            all_results[prop].append(entry)
    return all_results

def get_mean_std(all_results):
    all_results_mean_std = {}
    for prop in all_results:
        all_results_mean_std[prop] = []  # numpy array, size: (num_cycles, 3)
    for prop, result_all_cycles in all_results.items():
        for result_one_cycle in result_all_cycles:
            entry = [result_one_cycle[0], np.mean(result_one_cycle[-1]), np.std(result_one_cycle[-1])]
            all_results_mean_std[prop].append(entry)
        all_results_mean_std[prop] = np.array(all_results_mean_std[prop])
    return all_results_mean_std


def plot_results(all_results_mean_std):
    plt.figure(figsize=(10, 8))
    plt.text(0, 1.47, r'$S=8,  \, r=1$', fontsize=20)
    for prop, data in all_results_mean_std.items():
        plt.errorbar(x=data[:, 0], y=data[:, 1], yerr=data[:, 2], label=rf'$\sigma={prop} \rm \, \AA$',
                linewidth=2, capsize=6, elinewidth=2, markeredgewidth=2, fmt="-o", markersize=6)
        np.savetxt(f'saved/S8/sigma_{prop}.x', data[:, 0])
        np.savetxt(f'saved/S8/sigma_{prop}.y', data[:, 1])
        np.savetxt(f'saved/S8/sigma_{prop}.yerr', data[:, 2])
    plt.xlabel(r"$N \rm _{cycle}$", fontsize=20)
    plt.ylabel(r'$N \rm _{bond/atom}$', fontsize=20)
    plt.ylim(top=1.48)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18, frameon=False, loc='lower right')
    plt.show()

example_trj = f'analysis/S8/S_8_f_195_i_200_mi_1.0_s_13_sgm_3_cut_6_c_10_t_50000_at_2000/bond_per_atoms_selected.out'
num_cycles = list(map(int, pd.read_csv(example_trj, header=None)[0].to_numpy() / 200000))

num_incs = [3, 5, 7, 9, 11, 13]
seeds = [13, 14, 15]

all_results = get_results(num_cycles, num_incs, seeds)
all_results_mean_std = get_mean_std(all_results)
# print(all_results_mean_std)
plot_results(all_results_mean_std)