

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

# Local:
from functions import *


#========================================================================


# MAIN BODY

def run(args):

    PathToInputAln = args.i
    Partition = args.p
    Alpha = args.a

    if Partition:
        ListOfDicts = CodonSplitter(ReadSeq(PathToInputAln))
    else:
        ListOfDicts = ReadSeq(PathToInputAln)

    for DictName, SeqDict in enumerate(ListOfDicts):
        
        if Partition:
            print('\n')
            print(f"Performing all three tests on codon{DictName+1}")
        else:
            print(f"Performing all three tests on unpartitioned alignment")
            SeqDict = ListOfDicts

        AllPairs = list(itertools.combinations(SeqDict.keys(), 2))

        AllBowkers = []
        AllStuarts = []
        AllAbabnehs = []

        #paircount = 0
        nancount = 0
        for pair in tqdm(AllPairs):
            #paircount += 1
            #print(paircount)
            
            x, y = SeqDict[pair[0]], SeqDict[pair[1]]
            m = DivergenceMtx(x, y)
            
            if np.all(m[np.triu_indices(4, 1)] == 0):
                print("Encountered pair with all-zero off-diagonals")
                nancount +=1
                print(f"Invalid-pairs count: {nancount}")
                AllBowkers.append(np.nan)
                AllStuarts.append(np.nan)
                AllAbabnehs.append(np.nan)
                continue

            #print("="*30)
            #print(x)
            #print("-"*30)
            #print(y)
            
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

        print('\n')
        print("All Bowkers/Maximal symmetry stats:")
        print(AllBowkersMtx)
        print('\n')
        print("All Stuarts/Marginal symmetry stats:")
        print(AllStuartsMtx)
        print('\n')
        print("All Ababneh/Internal symmetry stats:")
        print(AllAbabnehsMtx)
        print('\n')
        # Print the three matrices to screen

        if Alpha == 0:
            BowkersAlpha = SequentialBonferroni(AllBowkers)
            StuartsAlpha = SequentialBonferroni(AllStuarts)
            AbabnehsAlpha = SequentialBonferroni(AllAbabnehs)
        else:
            BowkersAlpha = StuartsAlpha = AbabnehsAlpha = Alpha

        print(f"Printing Clustermaps for all three tests...")
        MaskedHeatmap(AllBowkersMtx, BowkersAlpha, f"{PathToInputAln}Codon{DictName+1}_Bowkers.png")
        MaskedHeatmap(AllStuartsMtx, StuartsAlpha, f"{PathToInputAln}Codon{DictName+1}_Stuarts.png")
        MaskedHeatmap(AllAbabnehsMtx, AbabnehsAlpha, f"{PathToInputAln}Codon{DictName+1}_Ababnehs.png")
        print('\n')

        if Partition:
            print(f"All three tests complete for partition {DictName+1} of alignment.")
            print('\n')
        else:
            print(f"All three tests complete for unpartitioned alignment.")
            print('\n')
        
        print(f"Three clustermaps of all pairwise scores have been written to {PathToInputAln}.")
        print('\n')
        print("="*79)
        # Save three heatmaps to where the input alignment is

        if not Partition:
            break
        # loop over single dict while keeping enumerate() for paritioning

#========================================================================


# WARNING FILTERS

simplefilter("ignore", ClusterWarning)
# Supresses warning re. X1 is too close to X1.T. Using sns.clustermap to make use of its ability to do map-permutations. It probably doesn't see that many symmetric matrices.
simplefilter("ignore", UserWarning)
# Supresses warnig re. max y and min y values being same for the data series.


#========================================================================


def main():
    
    parser = argparse.ArgumentParser(description='Use this to perform SRH tests on an alignment.')
    
    parser.add_argument("-input", help="Relative path of input alignment", required=True, dest="i", type=str)
    
    parser.add_argument("-partition", help="If True, SRHClusterMapper will partition input data into three codon positions and perform SRH tests on them separately. Defaults to False.", default=False, dest="p", type=bool)
    
    parser.add_argument("-alpha", help="Significance value. If given a custom/arbitrary value (e.g. 0.05), SRHClusterMapper will not perform Sequential Bonferroni correction. By default behaviour, Sequential Bonferroni correction will be performed to seek a significance value lower than 0.05. i.e. Leaving this option to default  will result in more sequences passing the symmetry tests.", default=0, dest="a", type=float)
    
    parser.set_defaults(func=run)
    
    args=parser.parse_args()
    
    args.func(args)


if __name__ == '__main__':
    main()
