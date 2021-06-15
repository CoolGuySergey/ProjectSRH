# Triggers the entire project

from functions import *
import argparse
import itertools

def run(args):

    PathToInputAln = args.i

    seqDict = readseq(PathToInputAln)
    allpairs = list(itertools.combinations(seqDict.keys(), 2))
    # n sequences will have n(n-1)/2 comparisons
    allBowkers = []
    allStuarts = []
    allAbabnehs = []

    for pair in allpairs:
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

    print (Broadcast2Matrix(allBowkers, seqDict))
    print (Broadcast2Matrix(allStuarts, seqDict))
    print (Broadcast2Matrix(allAbabnehs, seqDict))

def main():
    parser = argparse.ArgumentParser(description='Use this to do SRH tests on an alignmentl')                                                    
    parser.add_argument("-i", help="relative path of input", required=True, dest="i", type=str)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
