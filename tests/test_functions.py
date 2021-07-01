import functions

def test_ReadSeq():
    assert functions.ReadSeq("example.fa")

def test_CodonSplitter():
    ExampleDict = functions.ReadSeq("example.fa")
    assert functions.CodonSplitter(ExampleDict)
