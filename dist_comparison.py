import os
from collections import defaultdict, Counter

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from ete3 import Tree

COLOURS = reversed(['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628'])


def remove_certain_leaves(tr, to_remove=lambda node: False):
    """
    Removes all the branches leading to leaves identified positively by to_remove function.
    :param tr: the tree of interest (ete3 Tree)
    :param to_remove: a method to check is a leaf should be removed.
    :return: void, modifies the initial tree.
    """
    tips = [tip for tip in tr if to_remove(tip)]
    for tip in tips:
        while tip.is_leaf():
            if tip.is_root():
                return None
            parent = tip.up
            parent.remove_child(tip)
            tip = parent
    while len(tr.children) == 1:
        tr = tr.children[0]
        tr.dist = 0
        tr.up = None
    for parent in tr.traverse('postorder'):
        # If the parent node has only one child now, merge them.
        if len(parent.children) == 1:
            child = parent.children[0]
            child.dist += parent.dist
            grandparent = parent.up
            grandparent.remove_child(parent)
            grandparent.add_child(child)
    return tr


def get_node_name(n):
    return sorted([_.name for _ in n])


def read_tree(tree_path):
    for f in (3, 2, 5, 1, 0, 3, 4, 6, 7, 8, 9):
        try:
            return Tree(tree_path, format=f)
        except:
            continue
    raise ValueError('Could not read the tree {}. Is it a valid newick?'.format(tree_path))


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--trees', type=str, nargs='+')
    parser.add_argument('--consensus_tree', required=True, type=str)
    parser.add_argument('--age_pdf', required=True, type=str)
    parser.add_argument('--age_tab', required=True, type=str)
    params = parser.parse_args()

    dfs, acrs = [], []

    trees = [read_tree(t) for t in params.trees]
    names = [os.path.splitext(os.path.basename(tree))[0] for tree in params.trees]
    common_tips = set.intersection(*[{_.name for _ in t} for t in trees])
    trees = [remove_certain_leaves(t, lambda tip: tip.name not in common_tips) for t in trees]

    node_count = Counter()
    for t in trees:
        for n in t.traverse():
            n.name = 'ROOT' if n.is_root() else '_'.join(get_node_name(n))
            node_count[n.name] += 1

    branch2dist = defaultdict(list)
    for t in trees:
        for n in t.traverse('postorder'):
            if node_count[n.name] < len(trees):
                parent = n.up
                parent.remove_child(n)
                for child in list(n.children):
                    child.dist += n.dist
                    parent.add_child(child)
            else:
                branch2dist[n.name].append(n.dist)

    age_df = pd.DataFrame(columns=names)
    cons_tree = trees[0]
    for n in cons_tree.traverse('preorder'):
        n.dist = sum(branch2dist[n.name]) / len(branch2dist[n.name])
        age_df.loc[n.name, :] = branch2dist[n.name]
        if not n.is_root():
            age_df.loc[n.name, :] += age_df.loc[n.up.name, :]

    cons_tree.write(outfile=params.consensus_tree, format_root_node=True, format=3)

    for col in age_df.columns:
        age_df[col] = age_df[col].max() - age_df[col]

    age_df['avg'] = age_df.sum(axis=1) / age_df.shape[1]
    age_df.to_csv(params.age_tab, sep='\t', header=True, index_label='node')

    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    fig, axs = plt.subplots()
    for name, col, size, alpha in zip(names, COLOURS, range(90, 0, -14), range(65, 100, 5)):
        axs.scatter(age_df['avg'], age_df[name], s=size, facecolors='none', color=col, alpha=alpha / 100, linewidth=1, label=name)
    axs.legend(scatterpoints=1, fontsize=6)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    axs.set_xlabel("avg node age")
    axs.set_ylabel("node age")
    plt.tight_layout()
    plt.savefig(params.age_pdf, dpi=300, papertype='a5', orientation='landscape')
