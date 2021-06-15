from functions import *
import argparse
import itertools
from tqdm import tqdm

def run(args):

    PathToInputAln = args.i
    #for testing:
    #PathToInputAln = "MafftAligned/out.fa"
    
    seqDict = readseq(PathToInputAln)
    allpairs = list(itertools.combinations(seqDict.keys(), 2))
    # n sequences will have n(n-1)/2 comparisons
    allBowkers = []
    allStuarts = []
    allAbabnehs = []

    for pair in tqdm(allpairs):
        x,y = seqDict[pair[0]], seqDict[pair[1]]
        m = DivergenceMtx(x, y)
        BowkersStat, BowkersDf = list(Bowkers(m))
        BowkersPval = pval(BowkersStat, BowkersDf)
        StuartsStat, StuartsDf = Stuarts(m), 3
        StuartsPval = pval(StuartsStat, StuartsDf)
        AbabnehsStat, AbabnehsDf = Ababnehs(BowkersStat, BowkersDf, StuartsStat), BowkersDf-3
        AbabnehsPval = pval(AbabnehsStat, AbabnehsDf)
        #print (f"Current pair: {pair}")
        #print (f"p-values B:{BowkersPval} S:{StuartsPval} A:{AbabnehsPval}")
        allBowkers.append(BowkersPval)
        allStuarts.append(StuartsPval)
        allAbabnehs.append(AbabnehsPval)

    allBowkersMtx = Broadcast2Matrix(allBowkers, seqDict)
    allStuartsMtx = Broadcast2Matrix(allStuarts, seqDict)
    allAbabnehsMtx = Broadcast2Matrix(allAbabnehs, seqDict)
    
    print("All Bowkers stats/Maximal symmetry stats:")
    print (allBowkersMtx)
    print('\n\n')
    print("All Stuarts/Marginal symmetry stats:")
    print (allStuartsMtx)
    print('\n\n')
    print("All Ababneh/Internal symmetry stats:")
    print (allAbabnehsMtx)
    # Print the three matrices to screen

    MaskedHeatmap(allBowkersMtx, f"{PathToInputAln}_Bowkers")
    MaskedHeatmap(allStuartsMtx, f"{PathToInputAln}_Stuartss")
    MaskedHeatmap(allAbabnehsMtx, f"{PathToInputAln}_Ababnehs")
    print('\n\n')
    print(f"Three heatmaps have been written to wherever input alignment is")
    # Save three heatmaps to where the input alignment is

def main():
    parser = argparse.ArgumentParser(description='Use this to do SRH tests on an alignmentl')                                                    
    parser.add_argument("-i", help="relative path of input alignment", required=True, dest="i", type=str)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
