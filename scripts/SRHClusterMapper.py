

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


# WARNING FILTERS

simplefilter("ignore", ClusterWarning)
# Supresses warning re. X1 is too close to X1.T. Using sns.clustermap to make use of its ability to do map-permutations. It probably doesn't see that many symmetric matrices.
simplefilter("ignore", UserWarning)
# Supresses warnig re. max y and min y values being same for the data series.


#========================================================================


# MAIN BODY

def run(args):

    PathToInputAln = args.i
    Partition = args.p
    Alpha = args.a
    Benchmark = args.b

    if Partition:
        ListOfDicts = CodonSplitter(ReadSeq(PathToInputAln))
    else:
        ListOfDicts = ReadSeq(PathToInputAln)

    for DictName, SeqDict in enumerate(ListOfDicts):
        
        if Partition:
            TestName = f"codon{DictName+1}"
        else:
            TestName = "unpartitioned"
            SeqDict = ListOfDicts

        print('\n')
        print(f"Performing all three tests on {TestName} alignment.")

        AllPairs = list(itertools.combinations(SeqDict.keys(), 2))

        AllBowkers = []
        AllStuarts = []
        AllAbabnehs = []

        # Start looping through sequence pairs to buil DivergenceMtx
        # All three tests based on same DivergenceMtx
        nancount = 0
        for pair in tqdm(AllPairs):
            
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

            BowkersStat, BowkersDf = list(Bowkers(m))
            BowkersPval = pval(BowkersStat, BowkersDf)
            AllBowkers.append(BowkersPval)
            
            StuartsStat, StuartsDf = Stuarts(m), 3
            StuartsPval = pval(StuartsStat, StuartsDf)
            AllStuarts.append(StuartsPval)
            
            AbabnehsStat, AbabnehsDf = Ababnehs(BowkersStat, BowkersDf, StuartsStat), BowkersDf-3
            AbabnehsPval = pval(AbabnehsStat, AbabnehsDf)
            AllAbabnehs.append(AbabnehsPval)
        # End looping through sequence pairs

        # Print three matrices to screen
        AllBowkersMtx = Broadcast2Matrix(AllBowkers, SeqDict)
        print('\n')
        print("All Bowkers/Maximal symmetry stats:")
        print(AllBowkersMtx)

        AllStuartsMtx = Broadcast2Matrix(AllStuarts, SeqDict)
        print('\n')
        print("All Stuarts/Marginal symmetry stats:")
        print(AllStuartsMtx)
        print('\n')

        AllAbabnehsMtx = Broadcast2Matrix(AllAbabnehs, SeqDict)
        print("All Ababneh/Internal symmetry stats:")
        print(AllAbabnehsMtx)
        print('\n')

        if Alpha == 0:
            BowkersAlpha = SequentialBonferroni(AllBowkers)
            StuartsAlpha = SequentialBonferroni(AllStuarts)
            AbabnehsAlpha = SequentialBonferroni(AllAbabnehs)
        else:
            BowkersAlpha = StuartsAlpha = AbabnehsAlpha = Alpha

        print(f"Printing Clustermaps for all three tests...")
        print('\n')
        BowkersAllClusterDF = MaskedCluster(AllBowkersMtx, BowkersAlpha, f"{PathToInputAln}_{TestName}_Bowkers.jpg")
        StuartsAllClusterDF = MaskedCluster(AllStuartsMtx, StuartsAlpha, f"{PathToInputAln}_{TestName}_Stuarts.jpg")
        AbabnehsAllClusterDF = MaskedCluster(AllAbabnehsMtx, AbabnehsAlpha, f"{PathToInputAln}_{TestName}_Ababnehs.jpg")
        print(f"All three tests complete for {testname} alignment.")
        print('\n')
        print(f"Three clustermaps of all pairwise scores have been written to {PathToInputAln}.")
        print('\n')
        print("="*79)

        #Leftover = BowkersAllClusterDF
        ListOfStatNames = ["Bowkers", "Stuarts", "Ababnehs"]
        ListOfDFs = BowkersAllClusterDF, StuartsAllClusterDF, AbabnehsAllClusterDF
        StatNameDict = dict(zip(ListOfStatNames, ListOfDFs))
        
        for StatName, Leftover in StatNameDict.items():
            print('\n')
            print(f'Extracting {StatName} clusters:')
            print('\n')
            # Write clusters to OutDir, mkdir if it is yet to exist
            Outdir = f"{PathToInputAln}_{TestName}_{StatName}Clusters"
            if not os.path.exists(Outdir):
                os.mkdir(Outdir)
            
            ClusterNo = 0
            while len(Leftover) > 4:
                ClusterNo += 1
                NewCluster, Leftover = ExtractCluster(Leftover, Benchmark)
                WriteCluster(NewCluster, SeqDict, f"{Outdir}/Cluster{ClusterNo}_{len(NewCluster)}Seqs.fasta")
                print(f"Wrote cluster containing {len(NewCluster)} seqs.")
                print(f"There are {len(Leftover)} seqs left.")
                print('\n')
            else:
                if len(Leftover) > 0:
                    RemovedSeqs = Leftover.columns.tolist()
                    print(f"Removed Leftover Seqs: {RemovedSeqs}")
            print(f"{StatName} cluster extraction complete.")
            print('\n')
            print("="*79)
                    # Same behaviour as removing failing quartets
                    # In case of small blocks @ upper left corner
        
        if not Partition:
            break
        # loop over single dict while keeping enumerate() for paritioning


#========================================================================


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('on', 'yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('off', 'no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

    
def main():

    parser = argparse.ArgumentParser(
        description = 'Use this to perform SRH tests on an alignment.',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-i", "--input", help="Relative path of input alignment", required=True, dest="i", type=str)

    parser.add_argument("-p", "--partition", help="If true, SRHClusterMapper will partition input data into three codon positions and perform SRH tests on them separately. Defaults to False.", type=str2bool, default=False, dest="p")
    
    parser.add_argument("-a", "--alpha", help="Significance value. If given a custom/arbitrary value (e.g. 0.05), SRHClusterMapper will not perform Sequential Bonferroni correction. By default behaviour, Sequential Bonferroni correction will be performed to seek a significance value lower than 0.05. i.e. Leaving this option to default  will result in more sequences passing the symmetry tests.", default=0, dest="a", type=float)

    parser.add_argument("-b", "--benchmark", help="Benchmark / minimal purity of clusters in float representation. SRHClusterMapper will write out clusters where at least {benchmark*100} percent of all pairwise comparisons are passing pairs. Default off, that is, leaving this option to default will result in no clusters being written out.", default=0, dest="b", type=float)
    
    parser.set_defaults(func=run)
    
    args = parser.parse_args()
    
    args.func(args)

#========================================================================


if __name__ == '__main__':
    main()
