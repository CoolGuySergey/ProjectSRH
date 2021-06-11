# Ababneh’s test between 2 sequences

def Ababnehs(BowkersStat,BowkersDF, StuartsStat):
    '''
    MaxSym_int test. The remaining of general symmetry breakage is attributable to breakage of internal symmetry. Ababneh's = Bowker's - Stuarts. < 0.05, there are changes in the relative substitution rates between the root-to-tip axis of the tree or between lineages. Obtaining the data by chance under global homogeneity (assumption III) is unlikely.
    (n -1)(n-2)/2 degrees of freedom
    '''

    if isinstance(BowkersStat,float) and isinstance(StuartsStat, float)==True:
        if (BowkersDF > 3):
            s = BowkersStat - StuartsStat
        else:
            return np.nan
    else:
        return np.nan
    return float(s)
