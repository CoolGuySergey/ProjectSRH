from ete3 import Tree


def read_tree(nwk):
    for format in (3, 0, 1, 2):
        try:
            return Tree(nwk, format=format)
        except:
            continue


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--tree', required=True, type=str)
    parser.add_argument('--age', required=True, type=float)
    parser.add_argument('--out_tree', required=True, type=str)
    parser.add_argument('--dates', required=True, type=str)
    params = parser.parse_args()

    tree = read_tree(params.tree)

    with open(params.dates, 'w+') as f:
        f.write('name,date\n')
        f.write('root,0\n')
        for _ in tree:
            f.write('{},{}\n'.format(_.name, params.age))
    tree.add_child(name='root', dist=0)
    tree.write(outfile=params.out_tree)
