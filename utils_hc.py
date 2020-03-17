import numpy as np
from nh2lgf import newick_to_edgelist


class Node:
    def __init__(self, name, parent, id_node, mutation_id, loss=False, tot_mutations=0, gt_build=True):
        self.name = name
        self.id_node = id_node
        self.parent = parent
        self.children = []
        self.loss = loss
        self.mut_id = mutation_id
        self.tot_mutations = tot_mutations
        if parent:
            parent.children.append(self)

            if gt_build:
                gt_par_cp = parent.genotype_profile.copy()
                if self.loss:
                    gt_par_cp[mutation_id] -= 1
                else:
                    gt_par_cp[mutation_id] += 1

                self.genotype_profile = gt_par_cp
        else:
            self.genotype_profile = [0 for x in range(tot_mutations)]

    def is_ancestor_of(self, node):
        par = node.parent
        while par:
            if par == self:
                return True
            par = par.parent
        return False

    def print_node_dot(self):
        if not self.parent == None:
            print('\t"%s" -> "%s";' % (self.parent.id_node, self.id_node))
            if not '-' in self.name:
                print('\t"{0}" [label="{1}"];'.format(self.id_node, self.name))
            else:
                print('\t"{0}" [color=indianred1, style=filled, label="{1}"];'.format(
                    self.id_node, self.name[:-1]))

    def print_node_dot_file(self, fout):
        if not self.parent == None:
            fout.write('\t"%s" -> "%s";\n' %
                       (self.parent.id_node, self.id_node))
            if not '-' in self.name:
                fout.write('\t"{0}" [label="{1}"];\n'.format(
                    self.id_node, self.name))
            else:
                fout.write('\t"{0}" [color=indianred1, style=filled, label="{1}"];\n'.format(
                    self.id_node, self.name[:-1]))


def add_edge(start, end):
    start.children.append(end)
    end.parent = start


def __print_tree(node):
    if len(node.children) == 0:
        node.print_node_dot()
    else:
        node.print_node_dot()
        for child in node.children:
            __print_tree(child)


def __print_tree_file(node, fout):
    if len(node.children) == 0:
        node.print_node_dot_file(fout)
    else:
        node.print_node_dot_file(fout)
        for child in node.children:
            __print_tree_file(child, fout)


def print_dot_tree(node):
    print('digraph phylogeny {')
    print('\t"{0}" [label="{1}"];'.format(node.id_node, node.name))
    __print_tree(node)
    print('}')


def print_dot_tree_file(node, fout):
    fout.write('digraph phylogeny {\n')
    fout.write('\t"{0}" [label="{1}"];\n'.format(node.id_node, node.name))
    __print_tree_file(node, fout)
    fout.write('}\n')


def __copy_tree_rec(node, cp_parent, nid_dict):
    node_cp = Node(node.name, cp_parent, node.id_node,
                   node.mut_id, loss=node.loss)
    nid_dict[node_cp.id_node] = node_cp
    if len(node.children) == 0:
        return
    for child in node.children:
        __copy_tree_rec(child, node_cp, nid_dict)


def copy_tree(root):
    nid_nodes = {}
    cp_root = Node(root.name, root.parent, root.id_node,
                   root.mut_id, tot_mutations=root.tot_mutations)
    nid_nodes[root.id_node] = cp_root
    for child in root.children:
        __copy_tree_rec(child, cp_root, nid_nodes)
    return cp_root, nid_nodes

# True if col1 contains col2


def contains(col1, col2):
    for i in range(len(col1)):
        if not col1[i] >= col2[i]:
            return False
    return True


def build_tree_from_file(filepath, mutations_names, mutations_ids, tot_mutations):
    from subprocess import Popen, PIPE

    rb_tree = Popen(['ruby', 'tree', '-m', filepath], stdout=PIPE)
    stdout, _ = rb_tree.communicate()
    stdout = stdout.decode("utf-8").strip()
    # print(stdout)

    tree = stdout
    node_dict, edges = newick_to_edgelist(tree)

    building_dictionary = {}
    nid = 0

    # Get germline nid:
    for k in node_dict:
        if node_dict[k] == 'germline':
            nid = k

    root = Node('germline', None, nid, -1, tot_mutations=tot_mutations)
    building_dictionary[nid] = root

    for edge in edges:
        e, s = edge
        # print(s, e)

        x = None
        try:
            x = building_dictionary[s]
        except KeyError:
            x_column_index = int(node_dict[s][1:]) - 1
            if '---' in mutations_names[x_column_index]:
                loss = True
            else:
                loss = False
            x = Node(mutations_names[x_column_index],
                     None, s, mutations_ids[x_column_index], gt_build=False, loss=loss)
            building_dictionary[s] = x

        y = None
        try:
            y = building_dictionary[e]
        except KeyError:
            y_column_index = int(node_dict[e][1:]) - 1
            if '---' in mutations_names[y_column_index]:
                loss = True
            else:
                loss = False
            y = Node(mutations_names[y_column_index],
                     None, e, mutations_ids[y_column_index], gt_build=False, loss=loss)
            building_dictionary[e] = y
        # print(x,y)
        add_edge(x, y)

    # print_dot_tree(root)
    calculate_genotype_profile_subtree(root, building_dictionary)

    return (root, building_dictionary)


def calculate_genotype_profile_subtree(node, nid_dict):
    if node.parent:
        gt_par_cp = node.parent.genotype_profile.copy()
        if node.loss:
            gt_par_cp[node.mut_id] -= 1
        else:
            gt_par_cp[node.mut_id] += 1
        node.genotype_profile = gt_par_cp

    else:
        # This assumes that the root is correctly initiated at all 0s
        pass
    if len(node.children) == 0:
        return
    for child in node.children:
        calculate_genotype_profile_subtree(child, nid_dict)


def prune_and_reattach(node_prune, node_reattach, nid_dict):
    if node_prune.is_ancestor_of(node_reattach):
        return False
    node_prune.parent.children.remove(node_prune)
    node_prune.parent = node_reattach
    node_reattach.children.append(node_prune)

    # print('pr in')
    check_subtree_losses(node_reattach, nid_dict)

    # rebuild genotype profile of pruned subtree
    calculate_genotype_profile_subtree(node_reattach, nid_dict)
    # print('pr out')

    return True


def is_loss_valid(node, mut_id):
    par = node.parent
    while par:
        if par.mut_id == mut_id:
            return True
        par = par.parent
    return False


def is_already_lost(node, mut_id):
    par = node.parent
    while par:
        if par.loss == True and par.mut_id == mut_id:
            return True
        par = par.parent
    return False


def delete_node(node, nid_dict):
    parent = node.parent
    # node.parent = None
    parent.children.remove(node)
    for child in node.children:
        child.parent = parent
        parent.children.append(child)
    nid_dict.pop(node.id_node)
    # print('Deleted node: (%s, %d)' % (node.name, node.id_node))
    node = None


def check_subtree_losses(node, nid_dict):
    if node.loss:
        valid = is_loss_valid(node, node.mut_id)
        lost = is_already_lost(node, node.mut_id)

        if not valid or lost:

            # print('to delete', node.name, node.id_node)
            delete_node(node, nid_dict)

    if len(node.children) == 0:
        return
    for child in node.children:
        check_subtree_losses(child, nid_dict)


def import_ilp_out(filepath, k_dollo, mutation_names):
    in_matrix = np.genfromtxt(filepath, skip_header=0, delimiter=' ')
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

    imported_tree, imported_dict = build_tree_from_file(
        filepath, mut_names, mut_ids, len(mutation_names))

    return ((imported_tree, imported_dict), in_matrix)


def import_scs_input(filepath):
    mat = []
    with open(filepath, 'r') as fin:
        for line in fin:
            mat.append(
                [int(x) for x in line.strip().split()]
            )
    return mat
