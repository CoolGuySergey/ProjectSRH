import pytest
import functions
#========================================================================
# PREPROCESSING

def test_Preprocessing():
    
    assert list(functions.ReadSeq("Example.fa").keys()) == [
        '>Hycleus cichorii',
        '>Dermestes maculatus',
    ]

    with pytest.raises(IndexError):
        functions.ReadSeq("WrappedExample.fa")
    with pytest.raises(ValueError):
        functions.ReadSeq("AAExample.fa")

    assert len(functions.CodonSplitter(functions.ReadSeq("Example.fa"))) == 3
    
#========================================================================
# SYMMETRY TESTS


def test_SymmetryTests():
    
    with pytest.raises(ValueError):
        functions.DivergenceMtx("AACCGGTT", "AACCGGTT")
