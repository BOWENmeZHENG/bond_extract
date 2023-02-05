import utils
import numpy as np
import pandas as pd

example_trj = f'analysis/f_196_i_100_mi_1.0_s_12_sgm_5_cut_8_c_10_t_50000_at_2000/bond_per_atoms_selected.out'
num_cycles = list(map(int, pd.read_csv(example_trj, header=None)[0].to_numpy() / 200000))

num_incs = [100, 300, 500, 700, 900, 1100]
seeds = [12, 13, 14, 15, 16]

all_results = utils.get_results(num_cycles, num_incs, seeds)
all_results_mean_std = utils.get_mean_std(all_results)
# print(all_results_mean_std)
utils.plot_results(all_results_mean_std)