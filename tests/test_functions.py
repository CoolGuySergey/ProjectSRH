

# Description: Functions for generating SRH heatmaps for an alignment


#========================================================================


# IMPORTS

import pytest
from itertools import combinations
import sys
import os
sys.path.insert(1, '../scripts/')
from functions import *


#========================================================================


# PREPROCESSING

def test_Preprocessing():

    assert list(ReadSeq("Example.fa").keys()) == [
        '>Hycleus cichorii',
        '>Dermestes maculatus',
    ]

    assert list(ReadSeq("WrappedExample.fa").keys()) == [
        '>NC_039657',
        '>NC_037200',
    ]
    
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
    # Data from "Identifying Optimal Models of Evolution" Jermiin 2017
    
#========================================================================


# DATA VISUALISATION

def test_DataVisualisation():

    test = [0.13, 1e-05, 0.003, 0.05, 0.00099, 0.06, 0.035]
    assert SequentialBonferroni(test) == 0.003

    # For the list of test PVals, the first three are significant,
    # as PVal < Bonferroni corrected PVal
    
    # PVal   	test	inverse	      Bonferroni corrected PVal
    # 0.00001	1	7	      0.05/7 = 0.007143
    # 0.00099	2	6	      0.05/6 = 0.00833
    # 0.00300	3	5	      0.05/5 = 0.01
    # ---------------------------------------------------------
    # 0.03500	4	4	      0.05/4 = 0.0125
    # 0.05000	5	3	      0.05/3 = 0.016667
    # 0.06000	6	2	      0.05/2 = 0.025
    # 0.13000	7	1	      0.05/1 = 0.05

    ExDict = ReadSeq("RealExample.fa")
    AllPairs = list(combinations(ExDict.keys(), 2))
    AllBowkers = []
    for pair in AllPairs:
        x, y = ExDict[pair[0]], ExDict[pair[1]]
        m = DivergenceMtx(x, y)
        BowkersStat, BowkersDf = list(Bowkers(m))
        BowkersPval = pval(BowkersStat, BowkersDf)
        AllBowkers.append(BowkersPval)
    AllBowkersMtx = Broadcast2Matrix(AllBowkers, ExDict)
    assert len(AllBowkersMtx) == 4
