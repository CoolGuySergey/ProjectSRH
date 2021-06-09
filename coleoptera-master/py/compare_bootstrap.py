import logging

from ete3 import Tree

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--tree', required=True, type=str)
    parser.add_argument('--log', required=True, type=str)
    parser.add_argument('--tbe', required=True, type=str)
    parser.add_argument('--fbp', required=True, type=str)
    params = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    tbe = Tree(params.tbe, format=0)
    fbp = Tree(params.fbp, format=0)

    with open(params.log, 'w+') as f:
        f.write('node\tTBE\tFBP\tTBE-FBP\n')
        i = 0
        for n_tbe, n_fbp in zip(tbe.traverse('preorder'), fbp.traverse('preorder')):
            if n_tbe.is_leaf():
                continue
            n_fbp.support /= 100.0
            n_tbe.support = round(n_tbe.support, 2)
            f.write('{}\t{:.2f}\t{:.2f}\t{:.2f}\n'.format(i, n_tbe.support, n_fbp.support, n_tbe.support - n_fbp.support))
            n_tbe.support -= n_fbp.support
    tbe.write(outfile=params.tree, format=2)
