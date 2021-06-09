#Here I'm playing with bits of code to see what they do
import numpy as np
import math


def nCr(n,r):
    '''
    The factorial function
    '''
    f = math.factorial
    return f(n) // f(r) // f(n-r)  # nCr = n!/(r!(n-r)!)


def simMtx(x, y):
    '''
    divergence matrix
    '''
    a = 'ACGT'
    x = np.array(list(x))
    y = np.array(list(y))

    ax = (x[:, None] == a[None, :]).astype(int)
    ay = (y[:, None] == a[None, :]).astype(int)
    # array[:, None] smears array vertically
    # the other way round smears array horizontally
    # i.e. ax and ay will be 4 col wide and as tall as your seq is long
    # i.e. ax and ay  essentially represent the seq in a one-hot matrix
    
    return np.dot(ay.T, ax) # len*4.T into len*4 = 4*len into len*4
    # this will be you 4*4 divergence matrix, m

    
def MPTS(m):
    '''
    MaxSym test/Bowker's test. If < 0.05, general symmetry is violated. Obtaining the data by chance under stationarity (assumption I) or global homogeneity (assumption III) is unlikely.
    '''
    denominator = m+m.T
    # (dij + dji) all changes that are taking place
    # adding m to m.T makes denominator a 4*4  matrix
    # that's symmetrical across the main diagonal
    
    numerator = np.power(m-m.T,2)
    # (dij - dji)^2 squared difference between nucleotide i
    # becoming nucleotide j and nucleotide j becoming nucleotide i
    # multiplying m-m.T by itself makes numerator a 4*4 matrix
    # that's symmetrical across the main diagonal
    
    off_diag_indices = np.triu_indices(4,1)
    # np.triu_indices is tricky
    # but basically tiu stands for upper triangle and
    # if you have a 4*4 matrix a = np.arange(16).reshape(4, 4)
    # and set a[np.triu_indices(4)]=0
    # a will look like
    # array([[ 0,  0,  0,  0],
    #        [ 4,  0,  0,  0],
    #        [ 8,  9,  0,  0],
    #        [12, 13, 14,  0]])
    # here, set a[np.triu_indices(4,1)]=0 will turn a into
    # array([[ 0,  0,  0,  0],
    #        [ 4,  5,  0,  0],
    #        [ 8,  9, 10,  0],
    #        [12, 13, 14, 15]])
    # the 1 exludes the main diagonal from the valid indices
    # it makes sense to only take the upper off-diagonals

    numerator = np.squeeze(np.asarray(numerator[off_diag_indices]))
    denominator = np.squeeze(np.asarray(denominator[off_diag_indices]))
    # here we use off_diag_indices to actually take the numbers
    # np.asarray step is kinda redundant
    # np.squeeze returns the input array
    # but with all dimensions of length 1 removed
    # for example array([[[0],
    #                     [1],
    #                     [2]]])
    # will become array([0, 1, 2])

    if denominator[np.where(denominator != 0)[0]].size != 0:
        # look for indices for where denominator is not zero
        # and check that there actually ARE any non-zeros
        s = np.sum(np.divide(numerator[nonzeros], denominator[nonzeros]))
        # divide appropriate pairs and add them all up
        return float(s)
    else:
        return np.nan


def MPTS_df(m):
    '''
    degrees of freedom for the MaxSym test
    '''
    denominator = m+m.T
    off_diag_indices=np.triu_indices(len(denominator),1)
    denominator = np.squeeze(np.asarray(denominator[off_diag_indices]))
    i = np.count_nonzero(denominator)
    return int(i)
