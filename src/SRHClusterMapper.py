

# Description: Triggers entire project with functions contained in functions.py.


#========================================================================


# IMPORTS

# Standard:
import argparse
import itertools
from tqdm import tqdm

# Third party:
from scipy.cluster.hierarchy import ClusterWarning
from warnings import simplefilter
simplefilter("ignore", ClusterWarning)
# Supresses warning re. X1 is too close to X1.T. Using sns.clustermap to make use of its ability to do map-permutations. It probably doesn't see that many symmetric matrices.

# Local:
from functions import *


#========================================================================


# MAIN BODY

def run(args):

    PathToInputAln = args.i
    Partition = args.p

    if Partition:
        ListOfDicts = CodonSplitter(ReadSeq(PathToInputAln))
    else:
        ListOfDicts = ReadSeq(PathToInputAln)

    for DictName, SeqDict in enumerate(ListOfDicts):
        
        if Partition:
            print(f"Performing all three tests on codon{DictName+1}")
        else:
            print(f"Performing all three tests on unpartitioned alignment")
            SeqDict = ListOfDicts

        AllPairs = list(itertools.combinations(SeqDict.keys(), 2))

        AllBowkers = []
        AllStuarts = []
        AllAbabnehs = []

        for pair in tqdm(AllPairs):
            x, y = SeqDict[pair[0]], SeqDict[pair[1]]
            m = DivergenceMtx(x, y)
            BowkersStat, BowkersDf = list(Bowkers(m))
            BowkersPval = pval(BowkersStat, BowkersDf)
            StuartsStat, StuartsDf = Stuarts(m), 3
            StuartsPval = pval(StuartsStat, StuartsDf)
            AbabnehsStat, AbabnehsDf = Ababnehs(BowkersStat, BowkersDf, StuartsStat), BowkersDf-3
            AbabnehsPval = pval(AbabnehsStat, AbabnehsDf)
            AllBowkers.append(BowkersPval)
            AllStuarts.append(StuartsPval)
            AllAbabnehs.append(AbabnehsPval)

        AllBowkersMtx = Broadcast2Matrix(AllBowkers, SeqDict)
        AllStuartsMtx = Broadcast2Matrix(AllStuarts, SeqDict)
        AllAbabnehsMtx = Broadcast2Matrix(AllAbabnehs, SeqDict)

        print("All Bowkers stats/Maximal symmetry stats for:")
        print (AllBowkersMtx)
        print('\n\n')
        print("All Stuarts/Marginal symmetry stats:")
        print (AllStuartsMtx)
        print('\n\n')
        print("All Ababneh/Internal symmetry stats:")
        print (AllAbabnehsMtx)
        # Print the three matrices to screen

        print(f"Printing Clustermaps for all three tests...")
        MaskedHeatmap(AllBowkersMtx, f"{PathToInputAln}{DictName+1}_Bowkers.png")
        MaskedHeatmap(AllStuartsMtx, f"{PathToInputAln}{DictName+1}_Stuarts.png")
        MaskedHeatmap(AllAbabnehsMtx, f"{PathToInputAln}{DictName+1}_Ababnehs.png")
        print('\n')

        if Partition:
            print(f"All three tests complete for partition {DictName+1} of alignment.")
        else:
            print(f"All three tests complete for unpartitioned alignment.")
        
        print(f"Three clustermaps of all pairwise scores have been written to {PathToInputAln}.")
        print('\n')
        # Save three heatmaps to where the input alignment is


#========================================================================


def main():
    parser = argparse.ArgumentParser(description='Use this to do SRH tests on an alignment')                                                    
    parser.add_argument("-i", help="Relative path of input alignment", required=True, dest="i", type=str)
    parser.add_argument("-p", help="If True, SRHClusterMapper will partition input data into three codon positions and perform SRH tests on them separately. Defaults to False.", default=False, dest="p", type=bool)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
