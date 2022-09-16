#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import logging, warnings, h5py, multiprocess
import numpy as np
from cooler import balance, util, Cooler
from collections import defaultdict

# override the original functions

logger = logging.getLogger(__name__)

def _balance_genomewide(bias, cnv, ucnv, clr, spans, filters, map, tol, max_iters,
                        rescale_marginals, use_lock):
    scale = 1.0
    n_bins = len(bias)

    for _ in range(max_iters):
        marg = (
            balance.split(clr, spans=spans, map=map, use_lock=use_lock)
                .prepare(balance._init)
                .pipe(filters)
                .pipe(balance._timesouterproduct, bias)
                .pipe(balance._marginalize)
                .reduce(balance.add, np.zeros(n_bins))
        )
        
        nzmarg = marg[marg != 0]
        if not len(nzmarg):
            scale = np.nan
            bias[:] = np.nan
            var = 0.0
            break
        
        var_pool = []
        weights = []
        for uc in ucnv:
            if uc==0:
                continue
            umarg = marg[cnv==uc]
            tmp = umarg[umarg != 0]
            if not len(tmp):
                bias[cnv==uc] = 0
            else:
                umarg = umarg / tmp.mean()
                umarg[umarg==0] = 1
                bias[cnv==uc] = bias[cnv==uc] / umarg
                var_pool.append(tmp.var())
                weights.append(tmp.size)

        var = np.average(var_pool, weights=weights)
        logger.info("variance is {}".format(var))
        if var < tol:
            break
    else:
        warnings.warn(
            'Iteration limit reached without convergence.',
            balance.ConvergenceWarning)
    
    bias[bias==0] = np.nan
    for uc in ucnv:
        if uc==0:
            continue
        umarg = marg[cnv==uc]
        tmp = umarg[umarg != 0]
        if len(tmp):
            scale = tmp.mean()
            if rescale_marginals:
                bias[cnv==uc] = bias[cnv==uc] / np.sqrt(scale)

    '''
    # might be the reason why this normalized matrix has different range from ICE
    scale = nzmarg.mean()
    bias[bias == 0] = np.nan
    if rescale_marginals:
        bias /= np.sqrt(scale)
    '''
    return bias, var

def iterative_correction(clr, chunksize=int(1e7), map=map, tol=1e-5, min_nnz=10,
                         min_count=0, mad_max=5, ignore_diags=1, max_iters=200,
                         rescale_marginals=True, use_lock=False):
    
    # CNV must exists
    cnv = clr.bins()['CNV'][:].values
    # Divide the number of elements into non-overlapping chunks
    nnz = clr.info['nnz']
    if chunksize is None:
        chunksize = nnz
        spans = [(0, nnz)]
    else:
        edges = np.arange(0, nnz+chunksize, chunksize)
        spans = list(zip(edges[:-1], edges[1:]))

    # List of pre-marginalization data transformations
    base_filters = []
    if ignore_diags:
        base_filters.append(balance.partial(balance._zero_diags, ignore_diags))

    # Initialize the bias weights
    n_bins = clr.info['nbins']
    bias = np.ones(n_bins, dtype=float)

    # Drop bins with too few nonzeros from bias
    if min_nnz > 0:
        filters = [balance._binarize] + base_filters
        marg_nnz = (
            balance.split(clr, spans=spans, map=map, use_lock=use_lock)
                .prepare(balance._init)
                .pipe(filters)
                .pipe(balance._marginalize)
                .reduce(balance.add, np.zeros(n_bins))
        )
        bias[marg_nnz < min_nnz] = 0

    filters = base_filters
    marg = (
        balance.split(clr, spans=spans, map=map, use_lock=use_lock)
            .prepare(balance._init)
            .pipe(filters)
            .pipe(balance._marginalize)
            .reduce(balance.add, np.zeros(n_bins))
    )

    # Drop bins with too few total counts from bias
    if min_count:
        bias[marg < min_count] = 0

    # MAD-max filter on the marginals
    if mad_max > 0:
        offsets = clr._load_dset('indexes/chrom_offset')
        for lo, hi in zip(offsets[:-1], offsets[1:]):
            c_marg = marg[lo:hi]
            marg[lo:hi] /= np.median(c_marg[c_marg > 0])
        logNzMarg = np.log(marg[marg>0])
        med_logNzMarg = np.median(logNzMarg)
        dev_logNzMarg = balance.mad(logNzMarg)
        cutoff = np.exp(med_logNzMarg - mad_max * dev_logNzMarg)
        bias[marg < cutoff] = 0
    
    bias[cnv==0] = 0
    ucnv = np.unique(cnv)
    # Do balancing
    bias, var = _balance_genomewide(
            bias, cnv, ucnv, clr, spans, base_filters, map, tol, max_iters,
            rescale_marginals, use_lock)

    stats = {
        'tol': tol,
        'min_nnz': min_nnz,
        'min_count': min_count,
        'mad_max': mad_max,
        'ignore_diags': ignore_diags,
        'converged': var < tol,
        'var': var,
    }

    return bias, stats

def matrix_balance(cool_uri, nproc=1, chunksize=int(1e7), mad_max=5,
                min_nnz=10, min_count=0, ignore_diags=1, tol=1e-5,
                max_iters=1000):
    '''
    Perform separate matrix balancing for regions with different copy numbers
    and output the bias vector in the "sweight" column.
    '''
    cool_path, group_path = util.parse_cooler_uri(cool_uri)
    # Overwrite the existing sweight column
    with h5py.File(cool_path, 'r+') as h5:
        grp = h5[group_path]
        if 'sweight' in grp['bins']:
            del grp['bins']['sweight']
    
    clr = Cooler(cool_uri)
    
    try:
        if nproc > 1:
            pool = multiprocess.Pool(nproc)
            map_ = pool.imap_unordered
        else:
            map_ = map
        
        bias, stats = iterative_correction(
                clr,
                chunksize=chunksize,
                tol=tol,
                min_nnz=min_nnz,
                min_count=min_count,
                mad_max=mad_max,
                max_iters=max_iters,
                ignore_diags=ignore_diags,
                rescale_marginals=True,
                use_lock=False,
                map=map_)
    finally:
        if nproc > 1:
            pool.close()
    
    if not stats['converged']:
        logger.error('Iteration limit reached without convergence')
        logger.error('Storing final result. Check log to assess convergence.')
    
    with h5py.File(cool_path, 'r+') as h5:
        grp = h5[group_path]
        # add the bias column to the file
        h5opts = dict(compression='gzip', compression_opts=6)
        grp['bins'].create_dataset('sweight', data=bias, **h5opts)
        grp['bins']['sweight'].attrs.update(stats)
