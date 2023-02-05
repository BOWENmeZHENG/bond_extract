import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def parse(trj_name, line_limit=5000000):
    steps = []
    box_sizes = {}
    xyz = {}
    i = 0
    with open(trj_name) as f_trj:
        while i < line_limit:
            line = f_trj.readline()
            i += 1
            # Get step
            if line == 'ITEM: TIMESTEP\n':
                line = f_trj.readline()
                i += 1
                line_list = [int(num_str) for num_str in line.split()]
                step = line_list[0]
                steps.append(step)
                xyz[step] = []

            # Get atom number
            if line == 'ITEM: NUMBER OF ATOMS\n':
                line = f_trj.readline()
                i += 1
                line_list = [int(num_str) for num_str in line.split()]
                num_atoms = line_list[0]
            
            # Get box size
            if line == 'ITEM: BOX BOUNDS pp pp pp\n':
                line = f_trj.readline()
                i += 1
                line_list = [float(num_str) for num_str in line.split()]
                box_sizes[step] = line_list[-1] * 2

            # Get xyz
            if line[:11] == 'ITEM: ATOMS':
                line = f_trj.readline()
                i += 1
                line_list = [float(num_str) for num_str in line.split()]
                xyz[step].append(line_list[2: 5])
                for _ in range(num_atoms - 1):  
                    line = f_trj.readline()
                    i += 1         
                    line_list = [float(num_str) for num_str in line.split()]
                    if line_list[1] == 1.:
                        xyz[step].append(line_list[2: 5])
                xyz[step] = np.array(xyz[step])
            # End
            if line == "":
                break
    steps = np.array(steps)
    
    return steps, box_sizes, xyz


def comp_bond(frame, boxsize, bl):
    num_atoms = len(frame)
    bond_count = 0
    for i, atom1 in enumerate(frame):
        for atom2 in frame[i + 1:]:
            # filter out the impossible
            if boxsize - bl > abs(atom1[0] - atom2[0]) > bl or boxsize - bl > abs(atom1[1] - atom2[1]) > bl or boxsize - bl > abs(atom1[2] - atom2[2]) > bl:
                continue
            
            # solve the boundary cases
            if atom1[0] - atom2[0] > boxsize - bl:
                atom2[0] += boxsize
            if atom1[0] - atom2[0] < -(boxsize - bl):
                atom2[0] -= boxsize
            if atom1[1] - atom2[1] > boxsize - bl:
                atom2[1] += boxsize
            if atom1[1] - atom2[1] < -(boxsize - bl):
                atom2[1] -= boxsize
            if atom1[2] - atom2[2] > boxsize - bl:
                atom2[2] += boxsize
            if atom1[2] - atom2[2] < -(boxsize - bl):
                atom2[2] -= boxsize

            dist = comp_dist(atom1, atom2)
            if dist < bl:
                bond_count += 1
    bond_per_atom = bond_count / num_atoms
    return bond_count, num_atoms, bond_per_atom

def comp_dist(xyz_1, xyz_2):
    return np.linalg.norm(xyz_2 - xyz_1)

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
                trj_name = f'f_196_i_200_mi_1.0_s_{seed}_sgm_{prop}_cut_{prop+3.0}_c_10_t_50000_at_2000'
                data = pd.read_csv(f'analysis/{trj_name}/bond_per_atoms_selected.out', header=None)
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
    plt.figure(figsize=(8, 6))
    for prop, data in all_results_mean_std.items():
        plt.errorbar(x=data[:, 0], y=data[:, 1], yerr=data[:, 2], label=f'{prop} inclusions',
                linewidth=2, capsize=6, elinewidth=2, markeredgewidth=2, fmt="-o", markersize=6)
    plt.xlabel(r"$N \rm _{cycle}$", fontsize=20)
    plt.ylabel('bond per atom', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.show()