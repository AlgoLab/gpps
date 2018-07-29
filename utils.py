import numpy as np
from tree import *

def import_ilp_out(filepath, k_dollo, mutation_names):
    in_matrix = np.genfromtxt(filepath, skip_header=0, delimiter=' ')
    print(in_matrix.shape)
    mut_names = []
    mut_ids = []

    mut_index = 0
    for mut in mutation_names:
        mut_names.append(mut)
        mut_ids.append(mut_index)
        for i in range(k_dollo):
            mut_names.append('%s---' % mut)
            mut_ids.append(mut_index)
        mut_index += 1

    # print(mut_names)
    # print(mut_ids)

    imported_tree, imported_dict = build_tree_from_file(filepath, mut_names, mut_ids, len(mutation_names))

    # ---------- debug start
    # print_dot_tree(imported_tree)
    # print('\t'.join(mut_names))
    # for i in range(in_matrix.shape[0]):
    #     print('\t'.join(str(int(x)) for x in in_matrix[i,:]))

    # ----------- debug end
    return ((imported_tree, imported_dict), in_matrix)
    # return (build_tree_from_np(in_matrix, mut_names, mut_ids, len(mutation_names)), in_matrix)

def import_scs_input(filepath):
    mat = []
    with open(filepath, 'r') as fin:
        for line in fin:
            mat.append(
                [int(x) for x in line.strip().split()]
            )
    return mat

    