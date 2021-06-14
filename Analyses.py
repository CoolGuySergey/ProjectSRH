import itertools as ite
import numpy as np
from scipy.stats import chi2

def pval(s,df):
    '''
    p-value of the chi-square test
    '''
    if math.isnan(s)==False and df>0:
        p = 1.-float(chi2.cdf(s,df))
        return p
    else:
        return np.nan

# Don't forget to check if fasta is awk unwrapped to take up single lines
with open("MafftAligned/out.fa", "r") as filein:
    fasta = [i.split('\n') for i in filein.read().strip().split('\n\n')]

seqIDs = fasta[0][::2]
seqs = fasta[0][1::2]
seqDict = dict(zip(seqIDs, seqs))

def nCr(n,r):
    '''
    The factorial function
    '''
    f = math.factorial
    return f(n) // f(r) // f(n-r) # nCr = n!/(r!(n-r)!)

def Test_aln(aln,dset,dat):
    '''
    the matrix of all the pair-wise comparisons
    '''
    aln_array = np.array([list(rec) for rec in fasta[0]], dtype=object)
 
    no = nCr(len(aln),2)*3*len([len(v) for v in dat.charsets.keys()])+1
    
    p=np.empty([no,7],dtype='U22')
    p[0] = np.array(['dataset','Charset','Test','pvalue','d','Sp1','Sp2'])
    for n in tqdm(dat.charsets.keys()):
        for q in ite.combinations(list(range(len(aln))),2): #iterating over all taxa for sites
            m = simMtx('ACGT',aln_array[:,dat.charsets[n]][q[0]].tostring().upper().decode(),aln_array[:,dat.charsets[n]][q[1]].tostring().upper().decode())
            d = (np.sum(m)-sum(m[np.diag_indices(4)]))/np.sum(m)
            i = i+1

   # df for Stuarts test is 3
   # df for Ababneh's test is Bowkersdf-3
