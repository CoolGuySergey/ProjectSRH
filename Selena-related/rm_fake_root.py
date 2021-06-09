from ete3 import Tree


def read_tree(nwk):
    for fmt in (3, 0, 1, 2):
        try:
            return Tree(nwk, format=fmt)
        except:
            continue
    raise ValueError('Could not parse your tree {}'.format(nwk))


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_tree', required=True, type=str)
    parser.add_argument('--out_tree', required=True, type=str)
    parser.add_argument('--rootdate', required=True, type=str)
    parser.add_argument('--outgroup', required=True, type=str)
    params = parser.parse_args()

    tree = read_tree(params.in_tree)
    root_date = 0
    n = next(_ for _ in tree if _.name != 'root')
    while n:
        root_date -= n.dist
        n = n.up

    print("Initial root date is {}".format(root_date))

    for _ in tree:
        if _.name == 'root':
            parent = _.up
            parent.remove_child(_)
            if len(parent.children) == 1:
                child = parent.children[0]
                if parent.is_root():
                    tree = child
                    root_date += child.dist
                    tree.dist = 0
                    tree.up = None
                else:
                    grandparent = parent.up
                    grandparent.remove_child(parent)
                    grandparent.add_child(child, dist=parent.dist + child.dist)
            break

    print("Root date after removing the fake root is {}, tip number is {}".format(root_date, len(tree)))

    with open(params.outgroup, 'r') as f:
        outgroup = {_.strip() for _ in f.read().split('\n') if _.strip()}

    if outgroup:
        out_ns = [_ for _ in tree if _.name in outgroup]
        print('Outgroup is {}'.format(', '.join(_.name for _ in out_ns)))
        anc = out_ns[0] if len(out_ns) == 1 else out_ns[0].get_common_ancestor(out_ns)
        if len(anc) != len(out_ns):
            print('The outgroup is not monophyletic: {} tips vs {}:\n{}'.format(len(anc), len(out_ns), ', '.join(_.name for _ in anc)))
        if anc not in tree.children:
            print('There is potentially a rooting problem with the outgroup: it is not an immediate child of the root')
            while anc.up != tree:
                anc = anc.up
            print('Extended the outgroup to contain {}'.format(', '.join(_.name for _ in anc)))
        if len(tree.children) == 2:
            child = next(_ for _ in tree.children if _ != anc)
            root_date += child.dist
            tree = child
            tree.up = None
            tree.dist = 0
        else:
            tree.remove_child(anc)

        print("Root date after removing the outgroup is {}, tip number is {}".format(root_date, len(tree)))

    tree.write(outfile=params.out_tree)
    with open(params.rootdate, 'w+') as f:
        f.write('{}'.format(root_date))
