# Triggers the entire project

from functions import *
import argparse

#PathToInputAln = "MafftAligned/out.fa"

def run(args):

    PathToInputAln = args.i

    x = list(readseq(PathToInputAln).values())[0]
    y = list(readseq(PathToInputAln).values())[1]

    m = DivergenceMtx(x, y)

    BowkersStat, BowkersDf = list(Bowkers(m))
    BowkersPval = pval(BowkersStat, BowkersDf)

    StuartsStat, StuartsDf = Stuarts(m), 3
    
    StuartsPval = pval(StuartsStat, StuartsDf)
    AbabnehsStat, AbabnehsDf = Ababnehs(BowkersStat, BowkersDf, StuartsStat), BowkersDf-3
    AbabnehsPval = pval(AbabnehsStat, AbabnehsDf)

    print (BowkersPval, StuartsPval, AbabnehsPval)


def main():
    parser = argparse.ArgumentParser(description='Use this to do SRH tests on an alignmentl')                                                    
    parser.add_argument("-i", help="relative path of input", required=True, dest="i", type=str)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
