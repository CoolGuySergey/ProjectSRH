# Ababneh’s test between 2 sequences

# The remaining of general symmetry breakage is attributable to breakage of internal symmetry. Ababneh's = Bowker's - Stuarts

# (n -1)(n-2)/2 degrees of freedom, where n is the number of categories. If S_I^2 < 0.05, there are changes in the relative substitution rates between the root-to-tip axis of the tree or between lineages. Obtaining the data by chance under global homogeneity (assumption III) is unlikely.

def MPTIS(MPTSs,MPTSDF, MPTMSs):
    '''
    MaxSym_int test
    '''
    if isinstance(MPTSs,float) and isinstance(MPTMSs,float)==True:
        if (MPTSDF > 3):
            s = MPTSs-MPTMSs
        else:
            return np.nan
    else:
        return np.nan
    return float(s)
