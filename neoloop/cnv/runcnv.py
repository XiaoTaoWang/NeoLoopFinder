#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import os, cooler, pyBigWig
from cooler import balance
import numpy as np
import pandas as pd
from collections import Counter

def get_marginals(uri, exclude=['M', 'Y', 'MT', 'EBV'], chunksize=int(1e7), nproc=1):

    clr = cooler.Cooler(uri)

    if nproc > 1:
        pool = balance.Pool(nproc)
        map_ = pool.imap_unordered
    else:
        map_ = map

    nnz = clr.info['nnz']
    n_bins = clr.info['nbins']
    edges = np.arange(0, nnz+chunksize, chunksize)
    spans = list(zip(edges[:-1], edges[1:]))

    marg = (
        balance.split(clr, spans=spans, map=map_, use_lock=False)
            .prepare(balance._init)
            .pipe([])
            .pipe(balance._marginalize)
            .reduce(balance.add, np.zeros(n_bins))
    )
    table = clr.bins()[:][['chrom', 'start', 'end']]
    table['Coverage'] = marg.astype(int)
    pool = []
    chroms = [c for c in clr.chromnames if ((not c.lstrip('chr') in exclude) and (not '_' in c))]
    for chrom in chroms:
        pool.append(table[table['chrom']==chrom])
    
    table = pd.concat(pool)

    return table, clr.binsize

def signal_from_bigwig(table, bw_fil, name='GC'):
    
    bw = pyBigWig.open(bw_fil)
    arr = []
    for i in table.index:
        row = table.loc[i]
        v = bw.stats(row['chrom'], row['start'], row['end'], type='mean')[0]
        if v is None:
            v = 0
        arr.append(v)
    
    table[name] = np.r_[arr]

    return table

def count_REsites(table, npz_fil, res):

    RE = np.load(npz_fil)

    RE_by_bin = {}
    for chrom in RE:
        tmp = RE[chrom] // res
        RE_by_bin[chrom] = Counter(tmp)
    
    arr = []
    for i in table.index:
        row = table.loc[i]
        b_i = row['start'] // res
        chrom = row['chrom']
        arr.append(RE_by_bin[chrom][b_i])
    
    table['RE'] = np.r_[arr]

    return table

def filterZeros(table, cols=['GC', 'Mappability', 'RE']):

    mask = (table['GC'] != 0) & (table['Mappability'] != 0) & (table['RE'] != 0) & (table['Coverage'] != 0)
    filtered = table[mask]

    return mask, filtered


    
