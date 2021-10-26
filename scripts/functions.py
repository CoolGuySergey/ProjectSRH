	

# Description: Functions for generating SRH heatmaps for an alignment


#========================================================================


# IMPORTS

# Standard:
import itertools as ite
import math
import os
import sys

# Third party:
import numpy as np
from scipy.stats import chi2
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
matplotlib.use('Agg')

# Matplotlib is not thread-safe: in fact
# Threads must be set up the proper locks to serialize access to Matplotlib artists.
# Agg (non-interactive backend) so that one can work on separate figures from separate threads
# most GUI backends require being run from the main thread as well


#========================================================================


# PREPROCESSING
# Functions for reading in fasta (wrapped or unwrapped)
# Functions for partitioning input alignment


def ReadSeq(Path):

    """
    Reads in alignment as fasta file.
    
    In: (1 item) String where string contains relative path to alignment/fasta.
    Out: (1 item) Dictionary where keys are SeqIDs and items are Seqs.
    """
    
    with open(str(Path), "r") as FileIn:
        Fasta = [i.split('\n') for i in FileIn.read().strip().split('\n\n')][0]
    
    SeqDict = {}
    CurrentID = None
    CurrentSeq = ''
    
    for Line in Fasta:
        
        Line = Line.strip()
        
        if len(Line) == 0:
            continue
        if (Line[0] == '>'):
            CurrentID = Line
            CurrentSeq = ''
            SeqDict[CurrentID] = CurrentSeq
        else:
            CurrentSeq += Line.upper()

        SeqDict[CurrentID] = CurrentSeq.upper()
    
    return SeqDict


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
    
    if x == y:
        raise ValueError("Sorry. Caught duplicate sequences.")
    
    x = np.array(list(x))
    y = np.array(list(y))
    a = np.array(list('ACGT'))
    
    ax = (x[:, None] == a[None, :]).astype(int)
    ay = (y[:, None] == a[None, :]).astype(int)
    
    # array[:, None] smears array vertically
    # the other way round smears array horizontally
    # i.e. ax and ay will be 4 col wide and as tall as seq is long
    # i.e. ax and ay essentially represent seq as one-hot matrix

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

    In: (1 item) Divergence matrix, m.
    Out: (1 item) Stuarts Stats as float.
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

    In: (1 item) Divergence matrix, m.
    Out: (1 item) Ababnehs Stats as float.
    '''

    if isinstance(BowkersStat,float) and isinstance(StuartsStat, float):
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


def SequentialBonferroni(StatsList):

    """
    Seeks appropriate significance level from multiple p-values.
    
    In: (1 item) List of p-values.
    Out: (1 item) Float representing largest of the signficant p-values.
    """
    
    StatsList = sorted(StatsList)

    Rank = 0
    while StatsList[Rank] < (0.05 / (len(StatsList) - Rank)):
        Rank += 1

    return StatsList[Rank-1]


def Broadcast2Matrix(StatsList, SeqDict):
    
    '''
    Broadcast string of p-values to dataframe.
    
    In: (2 items) List of p-values, dictionary where keys are seqnames.
    Out: (1 item) Pandas dataframe where cells are p-values and row/cols 
    seqnames.
    '''
    
    n = len(SeqDict)
    mat = np.eye(n) # create n*n zero matrix
    
    iuu = np.triu_indices(n, 1) # 1 to exlude main diagonal
    mat[iuu] = StatsList # project string to upper triangular

    mat = mat + mat.T - np.diag(np.diag(mat))
    # Project to lower triangular by adding the transpose and
    # subtracting the diagonal. Will be faster albeit less readable.

    np.fill_diagonal(mat, np.nan)
    
    df = pd.DataFrame(mat, columns=SeqDict.keys())
    df.index = SeqDict.keys()
    
    return df

# for sequences a,b,c,d,e
# itertools/allpairs/allstats ALWAYS comes out in this order:
#         {"a" : [score with b, score with c, score with d, score with e],
#          "b" : [score with c, score with d, score with e],
#          "c" : [score with d, score with e],
#          "d" : [score with e]


def MaskedCluster(Dataframe, Alpha):
    
    '''
    Mask dataframe based on alpha and do row-col permutation.
    
    In: (2 items) Dataframe to visualise, significance level alpha.
    Out: (2 item) Boolean dataframe and cg post-clustering.
    '''

    # Initialise
    Boolean = Dataframe < Alpha
    cmap = sns.diverging_palette(240, 10, n=2)
    cg = sns.clustermap(Boolean, cmap=cmap, method='complete', metric='hamming', yticklabels=1, xticklabels=1)

    # Aesthetics and plotting
    # Font sizes:
    cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xmajorticklabels(), fontsize=1.75)
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize=1.75)
    # Hide unnecessaries:
    cg.ax_row_dendrogram.set_visible(False) # Hide 'trees'
    cg.ax_col_dendrogram.set_visible(False) # Hide 'trees'
    cg.cax.set_visible(False) # Hide colour bar
    
    # Return Reordered Boolean Dataframe
    RowReord = Boolean.iloc[cg.dendrogram_row.reordered_ind]
    FullReord = RowReord[[list(RowReord.columns)[x] for x in cg.dendrogram_row.reordered_ind]]

    return FullReord, cg


#========================================================================


# EXTRACTING CLUSTERS
# Functions for extracting out clusters
# Functions for writing cluster seqs to new fasta


def ExtractCluster(AllClusterDF, Benchmark):

    '''
    Write out bottom right cluster to CSV file.

    In: (2 items) Reordered dataframe, benchmark as float.
    Out: (2 item) Dataframe containing bottom right cluster.
                  Dataframe containign all remaining seqs.
    '''

    # AllClusterDF now looks like this:
    #        >Seq1  >Seq3  >Seq2  >Seq4
    # >Seq1  False  False   True   True
    # >Seq3  False  False   True   True
    # >Seq2   True   True  False  False
    # >Seq4   True   True  False  False

    # AllClusterDF.iloc[-2, -1:] is the second last row, last column
    # i.e. first meaningful score at bottom right of clustermark
    # AllClusterDF.iloc[-3, -2:] is the next row up, excluding main diag
    
    FailCount = sum(sum(AllClusterDF.iloc[-4:, -4:].to_numpy()))/2
    # FailCount initialised to that of bottom-right quartet
    AllComps = 6
    # There has been 6 comparisons thus far within bottom-right quartet
    Latch = -4
    
    # Seek passing starter quartet
    StarterLatch = 0
    while FailCount > round((1-Benchmark)*AllComps):
        StarterLatch += 1
        FailCount = sum(sum(AllClusterDF.iloc[-4-StarterLatch:-StarterLatch, -4-StarterLatch:-StarterLatch].to_numpy()))/2
        #print(f"{StarterLatch} dropped-seqs above failing quartet.")

    # If the previous loop has been triggered
    # Remove seqs that causes failing starter quartets
    if StarterLatch > 0:
        RemovedSeqs = AllClusterDF.iloc[-StarterLatch:, -StarterLatch:].columns.tolist()
        AllClusterDF = AllClusterDF.iloc[:-StarterLatch, :-StarterLatch]
        print("Starter quartet failed benchmark.")
        print(f"Shuffled to starter quartet that passes benchmark.")
        print(f"Removed intervening seqs: {RemovedSeqs}")
        
    else:
        print("Starter quartet passed benchmark.")
        
    # Expand passing starter quartet
    print("Expanding starter quartet into cluster.")
    while FailCount <= round((1-Benchmark)*AllComps):
        Latch -= 1

        # Break if (final) cluster reached end of map
        # Break if leftover cluster is smaller than a quartet
        # (rejection happens in SRHClusterMapper.py)
        if Latch == -len(AllClusterDF)-1 or len(AllClusterDF) < 4:
            break
        
        CurrentRow = AllClusterDF.iloc[Latch, Latch+1:]
        NewFails = sum(CurrentRow)
        FailCount += NewFails       # Update FailCount
        AllComps += len(CurrentRow) # Update Cluster size
        #print(f"{-Latch-4} rows above starter quartet.")
        #print(f"{FailCount} fails and {AllComps} comps thus far.")
        #print(f"Max fail is {round((1-Benchmark)*AllComps)}.")
    
    # Latest CurrentRow is the row that failied
    # Nip off cluster before starting next iteration.
    ClusterDF = AllClusterDF.iloc[Latch+1: , Latch+1:]
    RemainingDF = AllClusterDF.iloc[:Latch+1, :Latch+1]

    return ClusterDF, RemainingDF

          
def WriteCluster(ClusterDF, SeqDict, WritePath):

    """
    Reads wanted seqs from Cluster, writes appropriate seqs to WritePath.
    
    In: (3 items) Dataframe of seqs in cluster
                  SeqDict containing seqs
                  WritePath to new fasta
    Out: (1 item) Fasta file with selected entries.
    """

    WantedSeqNames = ClusterDF.columns.tolist()
    
    FastaOut = open(WritePath, "a")
    for WantedSeq in WantedSeqNames:
        FastaOut.write(WantedSeq + '\n')
        FastaOut.write(SeqDict[WantedSeq] + '\n')
    FastaOut.close()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For Use in Dev:

#ExDict = ReadSeq("94Seq_Supermatrix.fasta")
#AllPairs = list(ite.combinations(ExDict.keys(), 2))
#AllBowkers = []
#for pair in AllPairs:
#    x, y = ExDict[pair[0]], ExDict[pair[1]]
#    m = DivergenceMtx(x, y)
#    BowkersStat, BowkersDf = list(Bowkers(m))
#    BowkersPval = pval(BowkersStat, BowkersDf)
#    AllBowkers.append(BowkersPval)
#AllBowkersMtx = Broadcast2Matrix(AllBowkers, ExDict)
#BowkersAlpha = SequentialBonferroni(AllBowkers)
#
#df, cg = MaskedCluster(AllBowkersMtx, BowkersAlpha)
#ClusterDict = {">Galeruca daurica": 27}


def DemarcateCluster(cg, df, ClusterDict, Filename):

    '''
    Demarcate clusters based on info in ClusterDict and save image.
    
    In: (4 items) cg previously generated.
                  df (fullReord) previously generated
                  ClusterDict containing anchor (upperleft cell) and size.
                  Filename of png image.
    Out: (1 item) Saves image.
    '''
    
    ax = cg.ax_heatmap
    for AnchorSeq, ClusterSize in ClusterDict.items():
        # Identify upper left corner of box
        AnchorPos = df.index.to_list().index(AnchorSeq)
        ax.add_patch(mpatches.Rectangle((AnchorPos, AnchorPos), ClusterSize, ClusterSize, fill=False, edgecolor='yellow', lw=0.5))

    cg.savefig(Filename, format="jpg", dpi=450)
