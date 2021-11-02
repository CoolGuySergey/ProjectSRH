

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
    BenchmarkList = args.b

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
        print("All Bowkers/Maximal symmetry tests complete.")
        AllBowkersMtx.to_csv(f"{PathToInputAln}_{TestName}_AllBowkers.csv", index=False)
        print('\n')
        
        AllStuartsMtx = Broadcast2Matrix(AllStuarts, SeqDict)
        print("All Stuarts/Marginal symmetry tests complete.")
        AllStuartsMtx.to_csv(f"{PathToInputAln}_{TestName}_AllStuarts.csv", index=False)
        print('\n')

        AllAbabnehsMtx = Broadcast2Matrix(AllAbabnehs, SeqDict)
        print("All Ababneh/Internal symmetry tests complete.")
        AllAbabnehsMtx.to_csv(f"{PathToInputAln}_{TestName}_AllAbabnehs.csv", index=False)
        print('\n')

        if Alpha == 0:
            BowkersAlpha = SequentialBonferroni(AllBowkers)
            StuartsAlpha = SequentialBonferroni(AllStuarts)
            AbabnehsAlpha = SequentialBonferroni(AllAbabnehs)
        else:
            BowkersAlpha = StuartsAlpha = AbabnehsAlpha = Alpha

        BowkersAllClusterDF, BowkersCg = MaskedCluster(AllBowkersMtx, BowkersAlpha)
        StuartsAllClusterDF, StuartsCg = MaskedCluster(AllStuartsMtx, StuartsAlpha)
        AbabnehsAllClusterDF, AbabnehsCg = MaskedCluster(AllAbabnehsMtx, AbabnehsAlpha)

        print(f"All three tests complete for {TestName} alignment.")
        print('\n')
        print("="*79)

        ListOfStatNames = ["Bowkers", "Stuarts", "Ababnehs"]
        ListOfDFs = BowkersAllClusterDF, StuartsAllClusterDF, AbabnehsAllClusterDF
        ListofCgs = BowkersCg, StuartsCg, AbabnehsCg
        
        DFDict = dict(zip(ListOfStatNames, ListOfDFs))
        CgDict = dict(zip(ListOfStatNames, ListofCgs))
        
        for StatName, AllClusterDF in DFDict.items():            
            for Benchmark in BenchmarkList:
                # Reinitialise for next loop
                Leftover = AllClusterDF
                ClusterNo = 0
                ClusterDict = dict()
                CurImage = CgDict[StatName]
                CurDf = DFDict[StatName]
                
                print('\n')
                print(f'Extracting {StatName} clusters:')
                print(f'Current benchmark is {Benchmark}')
                print('\n')
                # Write clusters to OutDir, mkdir if it is yet to exist
                Outdir = f"{PathToInputAln}_{TestName}_{StatName}_Clusters_B{Benchmark}"
                if not os.path.exists(Outdir):
                    os.mkdir(Outdir)

                while len(Leftover) >= 4:
                    ClusterNo += 1
                    NewCluster, Leftover = ExtractCluster(Leftover, Benchmark)

                    if len(NewCluster) >= 4:
                        WriteCluster(NewCluster, SeqDict, f"{Outdir}/Cluster{ClusterNo}_{len(NewCluster)}Seqs.fasta")
                        print(f"Wrote cluster containing {len(NewCluster)} seqs.")
                        print(f"There are {len(Leftover)} seqs left.")
                        print('\n')

                        # Build ClusterDict to demarcate on map
                        AnchorSeq = NewCluster.columns.tolist()[0]
                        ClusterDict[AnchorSeq] = len(NewCluster)
                    else:
                        # Do not accept leftover clusters smaller than 4s
                        RemovedSeqs = NewCluster.columns.tolist()
                        print(f"Removed Leftover Seqs: {RemovedSeqs}")
                    
                else:
                    if len(Leftover) > 0:
                        RemovedSeqs = Leftover.columns.tolist()
                        print(f"Removed Leftover Seqs: {RemovedSeqs}")
                print(f"{StatName} cluster extraction complete.")
                print('\n')

                print(f"Printing {StatName} clustermap for {TestName} alignment.")
                print('\n')
            
                # Mark all clusters on cg generated from before
                DemarcateCluster(CurImage, CurDf, ClusterDict, f"{Outdir}.jpg")
                print(f"Three clustermaps of all pairwise scores have been written to location of {PathToInputAln}.")
                print('\n')
                print("="*79)
                # Same behaviour as removing failing quartets
                # In case of small blocks @ upper left corner

        # Logging
        with open(f'{PathToInputAln}_{TestName}.txt', 'w') as SummFile:
            
            SummFile.write('SUMMARY OF SRH ANALYSIS:\n')
            SummFile.write('\n')
            
            SummFile.write(f'Sequences count: {len(BowkersAllClusterDF)}\n')
            SummFile.write(f'Invalid-pairs count: {nancount}\n')
            NumbOfSeqs = len(BowkersAllClusterDF)
            allcomps = (((NumbOfSeqs ** 2) - NumbOfSeqs) * 0.5) - nancount
            SummFile.write(f'Total number of valid-pairs: {int(allcomps)}\n')
            SummFile.write('\n')

            FailingBTally = 0.5 * (BowkersAllClusterDF.sum().sum() - NumbOfSeqs)
            SummFile.write(f'Pairs that fail Bowkers: {int(FailingBTally)} ({FailingBTally/allcomps:.2f}%)\n')
            SummFile.write(f'Significance Level: {BowkersAlpha}\n')
            SummFile.write('\n')
            
            FailingSTally = 0.5 * (StuartsAllClusterDF.sum().sum() - NumbOfSeqs)
            SummFile.write(f'Pairs that fail Stuarts: {int(FailingSTally)} ({FailingSTally/allcomps:.2f}%)\n')
            SummFile.write(f'Significance Level: {StuartsAlpha}\n')
            SummFile.write('\n')
            
            FailingATally = 0.5 * (AbabnehsAllClusterDF.sum().sum() - NumbOfSeqs)
            SummFile.write(f'Pairs that fail Ababnehs: {int(FailingATally)} ({FailingATally/allcomps:.2f}%)\n')
            SummFile.write(f'Significance Level: {AbabnehsAlpha}\n')
            SummFile.write('\n')

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

    parser.add_argument("-b", "--benchmark", nargs='*', help="Benchmark / minimal purity of clusters in float representation. SRHClusterMapper will write out clusters where at least {benchmark*100} percent of all pairwise comparisons are passing pairs. Default off, that is, leaving this option to default will result in no clusters being written out.", default=0, dest="b", type=float)
    # This is the correct way to handle accepting multiple arguments.
    # '*' == 0 or more.
    # use: $ python whatever.py -b 1234 2345 3456 4567
    # [1234, 2345, 3456, 4567]
    
    parser.set_defaults(func=run)
    
    args = parser.parse_args()
    
    args.func(args)

#========================================================================


if __name__ == '__main__':
    main()
