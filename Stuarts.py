# Stuart's test between 2 sequences

import numpy as np
import math

def DivergenceMtx(x, y):
    '''
    divergence matrix
    '''
    a = np.array(list('acgt'))
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

    
def Stuarts(m):
    '''
    MaxSym_mar test/Stuartx's test for marginal symmetry. If < 0.05, marginal symmetry is violated. Obtaining the data by chance under stationarity (assumption I) is unlikely i.e. One of the four types of nucleotides is being substituted more than it is substituting other nucleotides.
    n-1 degrees of freedom, where n is the number of categories
    '''
    # Stuart's test statistics = u.T * V^-1 * u
    # u being the sequence pair’s vector of marginal differences
    # u.T = (d1• – d•1, d2• – d•2, d3• – d•3)
    # V being the variance-covariance matrix
    
    r = np.zeros((3)) # array([0., 0., 0.])
    r[0]=np.sum(m[0])
    r[1]=np.sum(m[1])
    r[2]=np.sum(m[2])
    c = [sum(row[i] for row in m) for i in range(4)]
    d = [r[0]-c[0],r[1]-c[1],r[2]-c[2]]
    

    ut = np.array([[d[0],d[1],d[2]]])
    u = ut.transpose()
    V = np.zeros((3,3))
    for (i,j) in ite.product(range(0,3),range(0,3)):
        if i==j:
            V[i,j]=r[i]+c[i]-2*m[i][i]
        elif i!=j:
            V[i,j]=-(m[i,j]+m[j,i])
        # Vij = dj• + d•j - 2dii  <-- if i=j
        # Vij = -(dij+dji)        <-- if i!=j

    if np.linalg.matrix_rank(V) != V.shape[0]:
        return np.nan
    else:
        Vi=np.linalg.inv(V)
        s = (ut.dot(Vi)).dot(u)[0][0]
        return float(s)

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
