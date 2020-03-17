#!/usr/bin/env python
import gurobipy as gp
from gurobipy import GRB
import sys
import os
import errno
import math
from datetime import datetime
import argparse
import itertools

from utils_ilp import read_matrix_tab, expand_name


#==================================================================#
#========================= PREREQUISITES ==========================#
#==================================================================#

#--------------------------Parse Arguments-------------------------#
parser = argparse.ArgumentParser(description='gpps- ILP', add_help=True)

parser.add_argument('-f', '--file', action='store', type=str, required=True,
                    help='path of the input file.')
parser.add_argument('-k', action='store', type=int, required=True,
                    help='k-value of the selected model. Eg: Dollo(k)')
parser.add_argument('-t', '--time', action='store', type=int, default=0,
                    help='maximum time allowed for the computation. Type 0 to not impose a limit (default).')
parser.add_argument('-o', '--outdir', action='store', type=str, required=True,
                    help='output directory.')

parser.add_argument('-d', '--maxdel', action='store', type=int,
                    default=-1, help='maximum number of deletion allowed')

parser.add_argument('-b', '--falsepositive', action='store', type=float, required=True,
                    help='set -b False positive probability.')
parser.add_argument('-a', '--falsenegative', action='store', type=float, required=True,
                    help='set -a False negative probability.')
parser.add_argument('--mps', action='store_true',
                    help='This will output the model in MPS format instead of running the solver')

args = parser.parse_args()


max_gains = 1
max_losses = int(args.k)
alpha = float(args.falsenegative)
beta = float(args.falsepositive)

#----------------------Initialize program----------------------#
input_matrix = read_matrix_tab(args.file)
matrix_name = os.path.basename(args.file).split('.')[0]

# Fixed parameters
#num_samples = len(input_matrix)
# num_clones = int(num_mutations * args.clones)
num_clones = len(input_matrix)
num_mutations = len(input_matrix[0])
num_columns = num_mutations * (1 + max_losses)
max_error = 1

#print('Num samples: %d' % num_samples)
print('Num mutations: %d' % num_mutations)
print('Num clones: %d' % num_clones)

try:
    os.makedirs(args.outdir)
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(args.outdir):
        pass
    else:
        raise

filename = os.path.splitext(os.path.basename(args.file))[0]
outfile = os.path.join(args.outdir, filename)

#==================================================================#
#========================== GUROBI MODEL ==========================#
#==================================================================#
model = gp.Model('Parsimony Phylogeny Model')
model.setParam('Threads', 4)
if args.time != 0:
    model.setParam('TimeLimit', args.time)

#---------------------------------------------------#
#------------------- VARIABLES ---------------------#
#---------------------------------------------------#

# -----------Variable Y and B---------------

lalpha = math.log(alpha)
lbeta = math.log(beta)
l_alpha = math.log(1 - alpha)
l_beta = math.log(1 - beta)
n_ones = 0
n_zeros = 0

print(lalpha, lbeta, l_alpha, l_beta)


print('Generating variables I, F, w, P.')
I = {}
F = {}
fminus = {}
P = {}
# False positive should be only a few
FP = {}

for row_index, row in enumerate(input_matrix):
    I[row_index] = {}
    F[row_index] = {}
    fminus[row_index] = {}
    P[row_index] = {}
    FP[row_index] = {}
    for col_index, cell in enumerate(row):
        # P
        names = expand_name(str(col_index), 1, max_losses)
        for name in names:
            P[row_index][name] = model.addVar(vtype=GRB.BINARY, obj=0,
                                              name='P{0}-{1}'.format(row_index, name))
        if cell < 2:
            # I
            I[row_index][col_index] = model.addVar(vtype=GRB.BINARY, obj=0,
                                                   name='I{0}-{1}'.format(row_index, col_index))
            model.update()
            model.addConstr(I[row_index][col_index] == cell,
                            "constr_I{0}-{1}".format(row_index, col_index))
            # F and fminus. fminus is equal to 1-F
            if cell == 0:
                F[row_index][col_index] = model.addVar(vtype=GRB.BINARY, obj=lalpha,
                                                       name='F{0}-{1}'.format(row_index, col_index))
                fminus[row_index][col_index] = model.addVar(vtype=GRB.BINARY, obj=l_beta,
                                                            name='f{0}-{1}'.format(row_index, col_index))
            if cell == 1:
                F[row_index][col_index] = model.addVar(vtype=GRB.BINARY, obj=l_alpha,
                                                       name='F{0}-{1}'.format(row_index, col_index))
                fminus[row_index][col_index] = model.addVar(vtype=GRB.BINARY, obj=lbeta,
                                                            name='f{0}-{1}'.format(row_index, col_index))
                FP[row_index][col_index] = model.addVar(vtype=GRB.BINARY, obj=0,
                                                        name='FP{0}-{1}'.format(row_index, col_index))
                model.update()
                model.addConstr(FP[row_index][col_index] == 1 - P[row_index][names[0]],
                                name='constr_FP{0}-{1}'.format(row_index, col_index))
            model.update()
            model.addConstr(F[row_index][col_index] == 1 - fminus[row_index][col_index],
                            name='constr_def_fminus{0}-{1}'.format(row_index, col_index))
            model.addConstr(F[row_index][col_index] == P[row_index][names[0]] - gp.quicksum(P[row_index][name] for name in names[1:]),
                            name='constr_balance_F{0}-{1}'.format(row_index, col_index))
            model.addConstr(P[row_index][names[0]] >= gp.quicksum(P[row_index][name] for name in names[1:]),
                            name='constr_imbalance_P{0}-{1}'.format(row_index, col_index))

# There are only a few false positives
model.addConstr(4 >= gp.quicksum(FP[r][c] for r in FP.keys() for c in FP[r].keys()),
                name='constr_few_FP')

model.update()

print('Generating variables B.')
B = {}
columns = list(itertools.chain.from_iterable(
    [expand_name(str(name), 1, max_losses) for name in range(num_mutations)]))
for p in columns:
    B[p] = {}
for p, q in itertools.combinations(columns, 2):
    B[p][q] = {}
    B[p][q]['01'] = model.addVar(vtype=GRB.BINARY, obj=0,
                                 name='B[{0},{1},0,1]'.format(p, q))
    B[p][q]['10'] = model.addVar(vtype=GRB.BINARY, obj=0,
                                 name='B[{0},{1},1,0]'.format(p, q))
    B[p][q]['11'] = model.addVar(vtype=GRB.BINARY, obj=0,
                                 name='B[{0},{1},1,1]'.format(p, q))
    model.update()
    for row_index, row in enumerate(input_matrix):
        model.addConstr(B[p][q]['01'] >= P[row_index][q] - P[row_index][p],
                        "constr_B01-{0}-{1}-{2}".format(p, q, row_index))
        model.addConstr(B[p][q]['10'] >= P[row_index][p] - P[row_index][q],
                        "constr_B10-{0}-{1}-{2}".format(p, q, row_index))
        model.addConstr(B[p][q]['11'] >= P[row_index][p] + P[row_index][q] - 1,
                        "constr_B11-{0}-{1}-{2}".format(p, q, row_index))
    model.addConstr(B[p][q]['01'] + B[p][q]['10'] + B[p][q]['11'] <= 2,
                    "constr_sum_B{0}-{1}".format(p, q))

deletions = {}
del_names = list()
for p in columns:
    if '-' in p:
        del_names.append(p)
        deletions[p] = model.addVar(
            vtype=GRB.BINARY, obj=0, name='Del[{0}]'.format(p))
        for row_index, row in enumerate(input_matrix):
            model.addConstr(deletions[p] >= P[row_index][p])

if args.maxdel != -1:
    model.addConstr(gp.quicksum(deletions[x]
                                for x in del_names) <= args.maxdel)

model.update()
model.modelSense = GRB.MAXIMIZE
model.update()

#---------------------------------------------------#
#-------------------- OPTIMIZE ---------------------#
#---------------------------------------------------#


if args.mps:
    print("Writing the model")
    model.write(outfile + '.mps')
    sys.exit(0)


print('#----- GUROBI OPTIMIZATION ----#')
start_optimize = datetime.now()
model.optimize()

#==================================================================#
#======================= POST OPTIMIZATION ========================#
#==================================================================#

matrix = []
# print(P)
with open('{0}.ilp.extended.out'.format(outfile), 'w+') as file_out:
    for row_index, row in enumerate(input_matrix):
        row_out = []
        for col_index, cell in enumerate(row):
            names = expand_name(str(col_index), 1, max_losses)
            for name in names:
                row_out.append(int(float(P[row_index][name].X)))
        matrix.append(row_out)
        file_out.write(' '.join([str(x) for x in row_out]))
        file_out.write('\n')

if model.status == GRB.Status.OPTIMAL or model.status == GRB.Status.TIME_LIMIT:
    value = float(model.objVal) + n_zeros * l_beta + n_ones * lbeta
    with open('{0}.ilp.log'.format(outfile), 'w+') as file_out:
        file_out.write('Optimal likelihood: %f\n' % value)
