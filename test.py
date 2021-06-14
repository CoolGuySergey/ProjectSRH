from functions import *

seqIDs = list("abcde")
seqsUpperCases = "accggatcgatcgagtcgagctga".upper(), "aggcttaggcggatgcgatcgtac".upper(), "accccggggcgggtttaaaccccg".upper(), "ggtcgagctaggcgatgcagcgat".upper(), "aggatcgagctgatcgatcgcccg".upper()
toyaln = dict(zip(seqIDs, seqsUpperCases))

import itertools

thing = list(itertools.combinations(toyaln.values(), 2))

[[(a[1]+b[1] if a[0]<b[0] else 0) for b in x] for a in x]

print("\t".join(['']+[a[0] for a in x]))
for a in x:
    print("\t".join([a[0]] + [(str(a[1]+b[1]) if a[0]<b[0] else '') for b in x]))
