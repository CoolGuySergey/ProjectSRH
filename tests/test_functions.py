import pytest
import functions
#========================================================================
# PREPROCESSING

def test_ReadSeq():
    assert functions.ReadSeq("Example.fa")
    with pytest.raises(FileNotFoundError):
        functions.ReadSeq("Null/or/Wrong/Path")
    with pytest.raises(AssertionError):
        functions.ReadSeq("AAExample.fa")

#========================================================================
# SYMMETRY TESTS
    
