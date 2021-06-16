from functions import *
import argparse
import itertools
from tqdm import tqdm
from scipy.cluster.hierarchy import ClusterWarning
from warnings import simplefilter
simplefilter("ignore", ClusterWarning)
# Supresses warning re. X1 is too close to X1.T. Using sns.clustermap to make use of its ability to do map-permutations. It probably doesn't see that many symmetric matrices

def run(args):

    PathToInputAln = args.i
    Partition = args.p
    
    #for testing:
    #PathToInputAln = "TransAligned/COX1.fa"
    #Partition = True
    #with open(PathToInputAln, "r") as filein:
    #    fasta = [i.split('\n') for i in filein.read().strip().split('\n\n')]
    #SeqIDs = fasta[0][::2]
    #SeqsOriginalCases = fasta[0][1::2]
    #SeqsUpperCases = [each_string.upper() for each_string in SeqsOriginalCases]

    if Partition == False:
        ListOfDicts = Readseq(PathToInputAln)
    else:
        ListOfDicts = CodonSplitter(Readseq(PathToInputAln))

    count=0
    for SeqDict in ListOfDicts:
        
        count +=1
        if Partition == False:
            print(f"Performing all three tests on unpartitioned alignment")
        else:
            print(f"Performing all three tests on codon{count}")

        AllPairs = list(itertools.combinations(SeqDict.keys(), 2))

        AllBowkers = []
        AllStuarts = []
        AllAbabnehs = []

        for pair in tqdm(AllPairs):
            x,y = SeqDict[pair[0]], SeqDict[pair[1]]
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

        MaskedHeatmap(AllBowkersMtx, f"{PathToInputAln}{count}_Bowkers.png")
        MaskedHeatmap(AllStuartsMtx, f"{PathToInputAln}{count}_Stuarts.png")
        MaskedHeatmap(AllAbabnehsMtx, f"{PathToInputAln}{count}_Ababnehs.png")
        print('\n')

        if Partition == False:
            print(f"All three tests complete for unpartitioned alignment.")
        else:
            print(f"All three tests complete for partition {count} of alignment.")
        print(f"Three clustermaps of all pairwise scores have been written to {PathToInputAln}.")
        print('\n')
        # Save three heatmaps to where the input alignment is

def main():
    parser = argparse.ArgumentParser(description='Use this to do SRH tests on an alignment')                                                    
    parser.add_argument("-i", help="Relative path of input alignment", required=True, dest="i", type=str)
    parser.add_argument("-p", help="If True, SRHClusterMapper will partition input data into three codon positions and perform SRH tests on them separately. Defaults to False.", default=False, dest="p", type=bool)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
