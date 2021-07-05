

# Description: Functions for generating SRH heatmaps for an alignment


#========================================================================


# IMPORTS

# Standard:
import itertools as ite
import math
import os

# Third party:
import numpy as np
from scipy.stats import chi2
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#========================================================================


# PREPROCESSING
# Functions for reading input alignment
# Functions for partitioning


def ReadSeq(path):

    """
    Reads in fasta file.
    
    In: (1 item) String where string contains relative path to fasta.
    Out: (1 item) Dictionary where keys are SeqIDs and items are Seqs.
    """

    with open(str(path), "r") as filein:
        fasta = [i.split('\n') for i in filein.read().strip().split('\n\n')]
    SeqIDs = fasta[0][::2]
    SeqsOriginalCases = fasta[0][1::2]
    SeqsUpperCases = [each_string.upper() for each_string in SeqsOriginalCases]

    if SeqIDs[1].startswith('>') == False:
        raise IndexError("""Sorry. Fasta sequence needs to be unwrapped first, try the following in your command line:
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < yourinput.fa > youroutput.fa
        """)
    if not set(SeqsUpperCases[0]).issubset(set("CGAT-")):
        raise ValueError('Sorry. Alignment seems to contain amino acids.')

    return dict(zip(SeqIDs, SeqsUpperCases))


def CodonSplitter(InputDict):
    
    """
    Partitions alignment into 1st/2nd/3rd codons.
    
    In: (1 item) Dictionary where keys are SeqIDs and items are Seqs.
    Out: (3 items) Dictionaries where keys are SeqIDs and items are Seqs
    identified as 1st/2nd/3rd codons.
    """
    
    PosOne = [each_string[::3] for each_string in InputDict.values()]
    PosTwo = [each_string[1::3] for each_string in InputDict.values()]
    PosThree = [each_string[2::3] for each_string in InputDict.values()]

    PosOneDict = dict(zip(InputDict.keys(), PosOne))
    PosTwoDict = dict(zip(InputDict.keys(), PosTwo))
    PosThreeDict = dict(zip(InputDict.keys(), PosThree))

    return PosOneDict, PosTwoDict, PosThreeDict


#========================================================================


# SYMMETRY TESTS
# Functions for build divergence matrix m between given seq-pair
# Funcitons for all three symmetry testsm, based on matrix m
# Functions for fetching p-values


def DivergenceMtx(x, y):

    '''
    Builds divergence matrix m between given seq-pair.
    
    In: (2 items) Strings where each string is seq in seq-pair.
    Out: (1 item) 4*4 numpy array representing divergence matrix, m.
    '''
    if x==y:
        raise ValueError("Sorry. Caught duplicate sequences.")
    
    x = np.array(list(x))
    y = np.array(list(y))
    a = np.array(list('ACGT'))
    
    ax = (x[:, None] == a[None, :]).astype(int)
    ay = (y[:, None] == a[None, :]).astype(int)
    
    # array[:, None] smears array vertically
    # the other way round smears array horizontally
    # i.e. ax and ay will be 4 col wide and as tall as your seq is long
    # i.e. ax and ay  essentially represent the seq in a one-hot matrix

    FinalMatrix = np.dot(ay.T, ax)
    # len*4.T into len*4 = 4*len into len*4
    
    return FinalMatrix


# First of three tests.
# df for Bowkers is returned in the same funciton
def Bowkers(m):
    
    '''
    MaxSym test/Bowker's test. If < 0.05, general symmetry is violated.
    Obtaining the data by chance under stationarity (assumption I) or
    global homogeneity (assumption III) is unlikely.
    
    n(n-1)/2 degrees of freedom, where n is the number of categories.
    
    In: (1 items) Divergence matrix, m.
    Out: (2 items) Bowkers Stats as float, degrees of freedom as integer.
    '''
    
    denominator = m + m.T
    # (dij + dji) all changes that are taking place
    # adding m to m.T makes denominator a 4*4  matrix
    # that's symmetrical across the main diagonal
    
    numerator = np.power(m - m.T, 2)
    # (dij - dji)^2 squared difference between nucleotide i
    # becoming nucleotide j and nucleotide j becoming nucleotide i
    # multiplying m-m.T by itself makes numerator a 4*4 matrix
    # that's symmetrical across the main diagonal
    
    OffDiagIndices = np.triu_indices(4, 1)
    # a = np.arange(16).reshape(4, 4)
    # a[np.triu_indices(4)]=0
    # a will look like
    # array([[ 0,  0,  0,  0],
    #        [ 4,  0,  0,  0],
    #        [ 8,  9,  0,  0],
    #        [12, 13, 14,  0]])
    #  a[np.triu_indices(4,1)]=0 will turn a into
    # array([[ 0,  0,  0,  0],
    #        [ 4,  5,  0,  0],
    #        [ 8,  9, 10,  0],
    #        [12, 13, 14, 15]])
    # the 1 exludes the main diagonal from the valid indices

    numerator = np.squeeze(np.asarray(numerator[OffDiagIndices]))
    denominator = np.squeeze(np.asarray(denominator[OffDiagIndices]))
    # here we use OffDiagIndices to actually take the numbers
    # np.asarray step is kinda redundant
    # np.squeeze returns the input array
    # but with all dimensions of length 1 removed
    # for example array([[[0],
    #                     [1],
    #                     [2]]])
    # will become array([0, 1, 2])
    
    nonzeros = np.where(denominator != 0)[0]
    #print(np.where(denominator != 0)[0])

    if np.count_nonzero(denominator) != 0:
        # check that I don't have an all-zero denominator
        s = np.sum(np.divide(numerator[nonzeros], denominator[nonzeros]))
        # divide appropriate pairs and add them all up
        df = np.count_nonzero(denominator)
        # how many types of changes are going on, should be 6 max
        # i.e. number of i,j pairs for which dij + dji > 0
        return float(s), int(df)
    else:
        return np.nan


# Second of three tests
def Stuarts(m):
    
    '''
    MaxSym_mar test/Stuartx's test for marginal symmetry. If < 0.05,
    marginal symmetry is violated. Obtaining the data by chance under
    stationarity (assumption I) is unlikely i.e. One of the four types of
    nucleotides is being substituted more than it is substituting other
    nucleotides.

    n-1 degrees of freedom, where n is the number of categories

    In: (1 items) Divergence matrix, m.
    Out: (2 items) Stuarts Stats as float.
    '''
    
    # Stuart's test statistics = u.T * V^-1 * u
    # u being the sequence pair’s vector of marginal differences
    # u.T = (d1• – d•1, d2• – d•2, d3• – d•3)
    # V being the variance-covariance matrix (3*3)
    
    r = np.zeros((3)) # array([0., 0., 0.])
    r[0] = np.sum(m[0])
    r[1] = np.sum(m[1])
    r[2] = np.sum(m[2])
    # summing across the rows of m, fetching dj•
    
    c = [sum(row[i] for row in m) for i in range(4)]
    # summing across the columns of m, fetching  d•j

    d = [r[0] - c[0], r[1] - c[1],r[2] - c[2]]
    
    ut = np.array([[d[0], d[1], d[2]]])
    u = ut.transpose()
    
    V = np.zeros((3, 3))
    for (i, j) in ite.product(range(0, 3),range(0, 3)):
        if i == j:
            V[i, j] = r[i] + c[i] - 2*m[i][i]
        elif i != j:
            V[i, j] = -(m[i,j] + m[j,i])
        # Vij = dj• + d•j - 2dii  <-- if i=j
        # Vij = -(dij+dji)        <-- if i!=j
        
    # V is now a 3*3 symmetric that generalises the notion of variance to multiple dimensions

    if np.linalg.matrix_rank(V) != V.shape[0]:
        # checking that there aren't any cols/rows that give us
        # no new information
        return np.nan
    else:
        Vi = np.linalg.inv(V)
        s = (ut.dot(Vi)).dot(u)[0][0] # u.T * V^-1 * u
        return float(s)


# Third of three tests
def Ababnehs(BowkersStat, BowkersDF, StuartsStat):
    
    '''
    MaxSym_int test. The remaining of general symmetry breakage is
    attributable to breakage of internal symmetry. Ababneh's = 
    Bowker's - Stuarts. < 0.05, there are changes in the relative
    substitution rates between the root-to-tip axis of the tree or
    between lineages. Obtaining the data by chance under global
    homogeneity (assumption III) is unlikely.
    
    (n-1)(n-2)/2 degrees of freedom. n can't be assumed to be 6, so this
    is calculated by subtracting 3 from BowkersDf.

    In: (1 items) Divergence matrix, m.
    Out: (2 items) Ababnehs Stats as float.
    '''

    if isinstance(BowkersStat,float) and isinstance(StuartsStat, float) == True:
        if (BowkersDF > 3):
            s = BowkersStat - StuartsStat
        else:
            return np.nan
    else:
        return np.nan
    return float(s)


# Fetch p-values
def pval(s, df):
    
    '''
    Fetches p-value of the chi-square test.

    In: (2 items) Test stat as float, degrees of freedom as integer.
    Out: (2 items) P-value as float.
    '''
    
    if math.isnan(s) == False and df > 0 :
        p = 1.-float(chi2.cdf(s, df))
        return p
    else:
        return np.nan


#========================================================================


# DATA VISUALISATION
# Functions for putting p-values into dataframe
# Functions for correcting significance threshold alpha
# Funcitons for building heatmap out of dataframe and threshold alpha


def Broadcast2Matrix(statsstring, seqDict):
    
    '''
    Broadcast string of p-values to dataframe.
    
    In: (2 items) String of p-values, dictionary where keys are seqnames.
    Out: (1 item) Pandas dataframe where cells are p-values and row/cols 
    seqnames.
    '''
    
    n = len(seqDict)
    mat = np.eye(n) # create n*n zero matrix
    iuu = np.triu_indices(n, 1) # 1 to exlude main diagonal
    mat[iuu] = statsstring # project string to upper triangular

    mat = mat + mat.T - np.diag(np.diag(mat))
    # Project to lower triangular by adding the transpose and
    # subtracting the diagonal. Will be faster albeit less readable.  
    
    df = pd.DataFrame(mat,columns=seqDict.keys())
    df.index = seqDict.keys()
    return df

# for sequences a,b,c,d,e
# itertools/allpairs/allstats ALWAYS comes out in this order:
#         {"a" : [score with b, score with c, score with d, score with e],
#          "b" : [score with c, score with d, score with e],
#          "c" : [score with d, score with e],
#          "d" : [score with e]


def MaskedHeatmap(dataframe, filename):
    
    '''
    Broadcast string of p-values to dataframe.
    
    In: (2 items) Dataframe to visualise, filename of png image.
    Out: (1 item) png image saved to working directory.
    '''
    
    alpha = 0.05 # uncorrected alpha level. No Bonferroni correction yet

    boolean = dataframe < alpha
    
    cmap = sns.diverging_palette(240,10,n=2)
    cg = sns.clustermap(boolean, cmap=cmap, yticklabels=1, xticklabels=1)
    cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xmajorticklabels(), fontsize=3)
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize=3)

    cg.ax_row_dendrogram.set_visible(False) # Hide 'trees'
    cg.ax_col_dendrogram.set_visible(False) # Hide 'trees'
    cg.cax.set_visible(False) # Hide colour bar
    cg.savefig(filename, format="png", dpi=250)
    #plt.show()
