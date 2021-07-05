

# Description: Functions for generating SRH heatmaps for an alignment


#========================================================================


# IMPORTS

import pytest
from functions import *


#========================================================================


# PREPROCESSING

def test_Preprocessing():
    
    assert list(ReadSeq("Example.fa").keys()) == [
        '>Hycleus cichorii',
        '>Dermestes maculatus',
    ]

    with pytest.raises(IndexError):
        ReadSeq("WrappedExample.fa")
    with pytest.raises(ValueError):
        ReadSeq("AAExample.fa")

    assert len(CodonSplitter(ReadSeq("Example.fa"))) == 3

    
#========================================================================


# SYMMETRY TESTS


def test_SymmetryTests():

    with pytest.raises(ValueError):
        DivergenceMtx("AACCGGTT", "AACCGGTT")

    ExDict = ReadSeq("RealExample.fa")
    ExDivergenceMtx = DivergenceMtx(ExDict[">Seq1"], ExDict[">Seq2"])
    BowkersStat, BowkersDf = list(Bowkers(ExDivergenceMtx))
    BowkersPval = pval(BowkersStat, BowkersDf)
    StuartsStat, StuartsDf = Stuarts(ExDivergenceMtx), 3
    StuartsPval = pval(StuartsStat, StuartsDf)
    AbabnehsStat, AbabnehsDf = Ababnehs(BowkersStat, BowkersDf, StuartsStat), BowkersDf-3
    AbabnehsPval = pval(AbabnehsStat, AbabnehsDf)
            
    assert BowkersPval == 2.5590610561176952e-05
    assert StuartsPval == 1.264401357392586e-06
    assert AbabnehsPval == 0.8500043787468997

    
#========================================================================


# DATA VISUALISATION
