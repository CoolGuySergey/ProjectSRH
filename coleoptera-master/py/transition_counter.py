from collections import Counter

import pandas as pd

from pastml.annotation import preannotate_forest
from pastml.tree import read_forest, annotate_dates
from pastml.acr import parse_date

STATE = 'pastml_state'
DATE = 'pastml_date'
LEVEL = 'pastml_level'


def state_combinations(node, columns, i=0):
    if i < len(columns):
        if i == len(columns) - 1:
            for state in getattr(node, columns[i]):
                yield (state,)
        else:
            for cmb in state_combinations(node, columns, i + 1):
                for state in getattr(node, columns[i]):
                    yield (state, *cmb)


def get_state_str(u):
    return u if len(u) > 1 else u[0]


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tree', help="named tree file produced by PastML, i.e. named.tree_....",
                        type=str, required=True)
    parser.add_argument('-d', '--root_date', required=False, default=None, type=str, nargs='*',
                        help="date(s) of the root(s) of the tree(s) (if not specified will be set to zero).")
    parser.add_argument('-a', '--acr', help="path to the ancestral character file produced by PastML, "
                                            "i.e. combined_ancestral_states.tab.", required=True, type=str)
    parser.add_argument('-c', '--columns', default=None, type=str, nargs='*',
                        help="columns of interest in the ancestral character file produced by PastML (by default all).")
    parser.add_argument('-o', '--output',
                        help="path to the output file that will contain transition dates.", required=True, type=str)
    parser.add_argument('-n', '--output_counts',
                        help="path to the output file that will contain transition counts.", required=False, type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.acr, header=0, index_col=0, sep='\t')
    df.index = df.index.map(str)
    forest = read_forest(params.tree)
    preannotate_forest(df, forest)

    if params.columns is None or len(params.columns) == 0:
        params.columns = df.columns

    for tree in forest:
        for n in tree.traverse():
            n.add_feature(STATE, set(state_combinations(n, params.columns)))

    root_dates = params.root_date
    if root_dates is not None:
        root_dates = [parse_date(d) for d in (root_dates if isinstance(root_dates, list) else [root_dates])]
        if 1 < len(root_dates) < len(forest):
            raise ValueError('{} trees are given, but only {} root dates.'.format(len(forest), len(root_dates)))
        elif 1 == len(root_dates):
            root_dates *= len(forest)
    annotate_dates(forest, date_feature=DATE, level_feature=LEVEL, root_dates=root_dates)

    data = []
    from_to_count = Counter()
    for tree in forest:
        for n in tree.traverse('preorder'):
            if n.is_root():
                continue
            states = getattr(n, STATE)
            up_states = getattr(n.up, STATE)
            for s in states:
                for u in up_states:
                    if s != u:
                        prob = 1 / len(states) / len(up_states)
                        data.append(
                            [get_state_str(u), get_state_str(s), getattr(n.up, DATE), getattr(n, DATE),
                             getattr(n, DATE) - getattr(n.up, DATE), prob, n.up.name, n.name])
                        from_to_count[(u, s)] += prob

    df = pd.DataFrame(data=data,
                      columns=['from state', 'to state', 'from date', 'to date', 'transition time', 'probability',
                               'from node id', 'to node id'])
    df.sort_values(by=['from date', 'transition time', 'from state', 'to state'], axis=0, inplace=True, ascending=False)
    df.to_csv(params.output, index=False, sep='\t')
    
    if params.output_counts:
        data = []
        for (u, s), count in from_to_count.items():
            data.append([get_state_str(u), get_state_str(s), count])
        df = pd.DataFrame(data=data, columns=['from state', 'to state', 'count'])
        df.sort_values(by=['count', 'from state', 'to state'], axis=0, inplace=True, ascending=False)
        df['percentage'] = 100 * df['count'] / df['count'].sum()
        df.to_csv(params.output_counts, index=False, sep='\t')
