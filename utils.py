import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def parse(trj_name, line_limit=10000000):
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

