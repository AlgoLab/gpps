# Copyright 2018
# Simomne Ciccolella
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from utils_hc import copy_tree, prune_and_reattach, build_tree_from_file, print_dot_tree_file, import_ilp_out, import_scs_input
from functools import lru_cache
import argparse
import os
import sys
import errno
import random


@lru_cache(maxsize=None)
def cell_row_likelihood(input_row, node_genotype, alpha, beta):
    likelihood = 0
    for j in range(len(input_row)):
        if input_row[j] == '0':
            if node_genotype[j] == '0':
                likelihood += np.log(1 - beta)
            elif node_genotype[j] == '1':
                likelihood += np.log(alpha)
            else:
                likelihood += -np.inf
        elif input_row[j] == '1':
            if node_genotype[j] == '0':
                likelihood += np.log(beta)
            elif node_genotype[j] == '1':
                likelihood += np.log(1 - alpha)
            else:
                likelihood += -np.inf

        elif input_row[j] == '2':
            likelihood += 0
    return likelihood


def greedy_tree_likelihood(tree, nid_dict, input_scs, alpha, beta):
    likelihood = 0
    attachment = []
    for row in input_scs:
        str_row = ''.join(str(int(x)) for x in row)
        best_lh = -np.inf
        best_attachment = -1
        for node in nid_dict:
            str_gt = ''.join(str(int(x))
                             for x in nid_dict[node].genotype_profile)

            lh = cell_row_likelihood(str_row, str_gt, alpha, beta)
            if lh > best_lh:
                best_lh = lh
                best_attachment = node
        likelihood += best_lh
        attachment.append(best_attachment)

    return likelihood, attachment


def get_expect_matrix(tree, nid_dict, input_scs, alpha, beta):
    _, attachment = greedy_tree_likelihood(
        tree, nid_dict, input_scs, alpha, beta)
    e_matrix = []

    for node in nid_dict:
        t_node = nid_dict[node]
        if t_node.loss and not t_node.id_node in attachment and len(t_node.children) == 0:
            print(t_node.name)

    for c in range(len(attachment)):
        e_matrix.append(nid_dict[attachment[c]].genotype_profile)
    return e_matrix


def generate_neighborhood(start_tree, start_nid_dict, neighborhood_size):
    neighbors = []
    while len(neighbors) < neighborhood_size:
        # prune-reattach only
        cp_tree, cp_dict = copy_tree(start_tree)
        node_ids = list(cp_dict.keys())
        prune = random.choice(node_ids)
        reattach = random.choice(node_ids)

        if prune != reattach:
            if prune_and_reattach(cp_dict[prune], cp_dict[reattach], cp_dict):
                neighbors.append((cp_tree, cp_dict))
            else:
                cp_tree = None
                cp_dict = None
    return(neighbors)


def hill_climbing(start_tree, start_nid_dict, neighborhood_size, max_iterations, alpha, beta, input_scs, proc_id=0):
    current_tree = start_tree
    current_dict = start_nid_dict
    # print('out_for')
    current_lh, _ = greedy_tree_likelihood(
        start_tree, start_nid_dict, input_scs, alpha, beta)
    print('start lh: %f' % current_lh)

    # print(input_scs)

    current_iteration = 0
    while current_iteration < max_iterations:

        print('Current iteration: %d' % current_iteration)

        neighbors = generate_neighborhood(
            current_tree, current_dict, neighborhood_size)
        next_eval = -np.inf
        next_sol = None

        for ng in neighbors:
            # print(ng)
            ng_lh, _ = greedy_tree_likelihood(
                ng[0], ng[1], input_scs, alpha, beta)

            if ng_lh > next_eval:
                next_eval = ng_lh
                next_sol = (ng[0], ng[1])

        if next_eval > current_lh:
            current_tree = next_sol[0]
            current_dict = next_sol[1]
            current_lh = next_eval
            print('Found a better solution with likelihood: %f' % current_lh)
        current_iteration += 1

    return current_tree, current_dict


parser = argparse.ArgumentParser(
    description='gpps- hill climber', add_help=True)
parser.add_argument('-i', '--ilpfile', action='store', type=str, required=True,
                    help='path of the ILP output file.')
parser.add_argument('-s', '--scsfile', action='store', type=str, required=True,
                    help='path of the SCS input file. (same input feeded to the ILP)')
parser.add_argument('-k', action='store', type=int, required=True,
                    help='k-value of the selected model. Eg: Dollo(k)')
parser.add_argument('-o', '--outdir', action='store', type=str, required=True,
                    help='output directory.')

parser.add_argument('-b', '--falsepositive', action='store', type=float, required=True,
                    help='set -b False positive probability.')
parser.add_argument('-a', '--falsenegative', action='store', type=float, required=True,
                    help='set -a False negative probability.')
parser.add_argument('--ns', action='store', type=int, required=True,
                    help='Hill climbing neighbourhood size.')
parser.add_argument('--mi', action='store', type=int, required=True,
                    help='Hill climbing maximum iterations.')
parser.add_argument('--names', action='store', type=str,
                    help='Mutation names.')
args = parser.parse_args()


FILEPATH_ILP = args.ilpfile
FILEPATH_SCS = args.scsfile
K_DOLLO = args.k
INPUT_SCS_MATRIX = import_scs_input(FILEPATH_SCS)
if not args.names:
    NAMES = [str(x) for x in range(1, len(INPUT_SCS_MATRIX[0]) + 1)]
else:
    with open(args.names, 'r') as fin:
        NAMES = []
        for line in fin:
            NAMES.append(line.strip())
OUT_DIR = args.outdir

OUTFILE = os.path.splitext(os.path.basename(args.ilpfile))[0]

built_tree, np_input = import_ilp_out(FILEPATH_ILP, K_DOLLO, NAMES)

ALPHA = args.falsenegative
BETA = args.falsepositive

try:
    os.makedirs(OUT_DIR)
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(OUT_DIR):
        pass
    else:
        raise

imported_tree, imported_nid_dict = built_tree
greedy_tree_likelihood(imported_tree, imported_nid_dict,
                       INPUT_SCS_MATRIX, ALPHA, BETA)

with open('%s/%s.input.gv' % (args.outdir, OUTFILE), 'w+') as fout:
    print_dot_tree_file(imported_tree, fout)
hc_best_tree, hc_best_dict = hill_climbing(imported_tree, imported_nid_dict, neighborhood_size=args.ns,
                                           max_iterations=args.mi, alpha=ALPHA, beta=BETA, input_scs=INPUT_SCS_MATRIX)


with open('%s/%s.hill_climbing.gv' % (args.outdir, OUTFILE), 'w+') as fout:
    print_dot_tree_file(hc_best_tree, fout)

with open('%s/%s.hill_climbing.scs.out' % (args.outdir, OUTFILE), 'w+') as fout:
    e_mat = get_expect_matrix(
        hc_best_tree, hc_best_dict, INPUT_SCS_MATRIX, ALPHA, BETA)
    for row in e_mat:
        fout.write(' '.join(str(x) for x in row))
        fout.write('\n')

# print(cell_row_likelihood.cache_info())
