#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import logging, copy
from collections import defaultdict
import numpy as np
from sklearn.mixture import GaussianMixture
from neoloop.cnv.loadcnv import binCNV
from pomegranate import NormalDistribution, HiddenMarkovModel

logger = logging.getLogger(__name__)

def cbs_stat(x):
    '''
    Given x, Compute the subinterval x[i0:i1] with the maximal segmentation statistic t. 
    Returns t, i0, i1.

    '''
    
    x0 = x - np.mean(x)
    n = len(x0)
    y = np.cumsum(x0)
    e0, e1 = np.argmin(y), np.argmax(y)
    i0, i1 = min(e0, e1), max(e0, e1)
    s0, s1 = y[i0], y[i1]

    return (s1-s0)**2*n/(i1-i0+1)/(n+1-i1+i0), i0, i1+1

def cbs(x, shuffles=1000, p=.05):
    '''
    Given x, find the interval x[i0:i1] with maximal segmentation statistic t. Test that statistic against
    given (shuffles) number of random permutations with significance p.  Return True/False, t, i0, i1; True if
    interval is significant, false otherwise.
    '''
    max_t, max_start, max_end = cbs_stat(x)
    if max_end-max_start == len(x):
        return False, max_t, max_start, max_end
    if max_start < 5:
        max_start = 0
    if len(x)-max_end < 5:
        max_end = len(x)
    thresh_count = 0
    alpha = shuffles*p
    xt = x.copy()
    for i in range(shuffles):
        np.random.shuffle(xt)
        threshold, s0, e0 = cbs_stat(xt)
        if threshold >= max_t:
            thresh_count += 1
        if thresh_count > alpha:
            return False, max_t, max_start, max_end

    return True, max_t, max_start, max_end

def rsegment(x, start, end, L=[], shuffles=1000, p=.05):
    '''
    Recursively segment the interval x[start:end] returning a list
    L of pairs (i,j) where each (i,j) is a significant segment.
    '''
    threshold, t, s, e = cbs(x[start:end], shuffles=shuffles, p=p)
    if (not threshold) | (e-s < 5) | (e-s == end-start):
        L.append((start, end))
    else:
        if s > 0:
            rsegment(x, start, start+s, L)
        if e-s > 0:
            rsegment(x, start+s, start+e, L)
        if start+e < end:
            rsegment(x, start+e, end, L)

    return L

def segment(x, shuffles=1000, p=1e-5):
    '''
    Segment the array x, using significance test based on shuffles rearrangements and significance level p
    '''
    start = 0
    end = len(x)
    L = []
    rsegment(x, start, end, L, shuffles=shuffles, p=p)

    return L

def combine_segment(sig, cbs_seg, hmm_seg, min_seg, min_diff, bufsize=4):

    start = cbs_seg[0][0]
    pool1 = set()
    for b in cbs_seg:
        for i in range(bufsize):
            pool1.add(b[1]-i)
            pool1.add(b[1]+i)

    pool2 = set([b[1] for b in hmm_seg])
    pool = [start] + sorted(pool1 & pool2)
    seg = []
    # combine segments with similar copy numbers
    start, end = pool[0], pool[1]
    if len(pool) > 2:
        for i in range(2, len(pool)):
            ext = pool[i]
            arr1 = sig[start:end]
            mask1 = arr1 == 0
            zero_ratio1 = mask1.sum() / mask1.size
            arr2 = sig[end:ext]
            mask2 = arr2 == 0
            zero_ratio2 = mask2.sum() / mask2.size
            if (mask1.size < min_seg) or (mask2.size < min_seg):
                _add_seg(seg, [start, ext])
                start, end = start, ext
            else:
                if (zero_ratio1 > 0.8) or (zero_ratio2 > 0.8):
                    _add_seg(seg, [start, end])
                    start, end = end, ext
                else:
                    v1 = np.log2(np.median(arr1[arr1!=0]))
                    v2 = np.log2(np.median(arr2[arr2!=0]))
                    if abs(v1 - v2) < min_diff:
                        _add_seg(seg, [start, ext])
                        start, end = start, ext
                    else:
                        _add_seg(seg, [start, end])
                        start, end = end, ext
                        
        if end > seg[-1][1]:
            _add_seg(seg, [start, end])
    else:
        _add_seg(seg, [start, end])

    return seg

def _add_seg(seg, r):

    if not len(seg):
        seg.append(r)
    else:
        if r[0]==seg[-1][0]:
            seg[-1] = r
        else:
            seg.append(r)

class HMMsegment(binCNV):

    def __init__(self, bedgraph, res, nproc=1, ploidy=2, n_states=None):

        binCNV.__init__(self, bedgraph, res)
        self.n_jobs = nproc
        self.ploidy = ploidy
        self.n_states = n_states

    def segment(self, min_seg=5, min_diff=0.4, p=1e-5, max_dist=4):

        outlines = []
        counts = defaultdict(int)
        vMap = {}
        training_seqs, gmm = self.get_states(maxdist=10)
        if self.n_states is None:
            n_states = len(gmm.means_.ravel())
            #logger.info('Estimated HMM state number: {0}'.format(n_states))
            #logger.info('Means of the hidden distributions: {0}'.format(gmm.means_.ravel()))
            #logger.info('Covariances of the hidden distributions: {0}'.format(gmm.covariances_.ravel()))
        else:
            n_states = self.n_states

        if n_states > 1:
            logger.info('Training HMM with {0} states...'.format(n_states))
            model = HiddenMarkovModel.from_samples(
                NormalDistribution,
                n_components=n_states,
                init='kmeans++',
                n_init=10,
                X=training_seqs,
                algorithm='baum-welch',
                stop_threshold=1e-3,
                max_iterations=1000,
                n_jobs=self.n_jobs
            )
            logger.info('Done')
        else:
            model = None

        for c in self.bin_cnv:
            logger.info('Segmenting Chromosome {0} ...'.format(c))
            sig = self.bin_cnv[c]
            if not model is None:
                queue = sig[sig > 0]
                seq = np.r_[[s.name for i, s in model.viterbi(np.log2(queue))[1][1:]]]
                hmm_seg = self.assign_cnv(queue, seq)
                predicted = np.zeros(sig.size)
                predicted[sig > 0] = hmm_seg
                hmm_seg = self.call_intervals(predicted)
            else:
                hmm_seg = [(0, sig.size)]

            #logger.info('Combine results from the CBS algorithm ...')
            logsig = np.zeros_like(sig)
            logsig[sig>0] = np.log2(sig[sig>0])
            cbs_seg = segment(logsig, p=p)
            seg = combine_segment(sig, cbs_seg, hmm_seg, min_seg,
                                  min_diff, bufsize=max_dist) # combine results from cbs and hmm
            
            for s, e in seg:
                tmp = sig[s:e]
                mask = tmp==0
                zero_ratio = mask.sum() / mask.size
                if zero_ratio > 0.8:
                    cn = 0
                else:
                    cn = int(np.round(np.median(tmp[tmp!=0])*self.ploidy))
                cur = [c, s, e, cn]
                outlines.append(cur)
                vMap[cn] = cn    
                counts[cn] += (cur[2]-cur[1])
        
        cn_min_count = min(counts.values())
        while cn_min_count < 5: # total number of bins with the a specific copy number must greater than this value
            for i in counts:
                if counts[i]==cn_min_count:
                    cn = i
            neighbor = 0
            mindis = 10000 # a very large number to start
            for i in counts:
                if i == cn:
                    continue
                if np.abs(i - cn) < mindis:
                    neighbor = vMap[i]
                    mindis = np.abs(i - cn)

            tc = counts[cn] + counts[neighbor]

            for i in vMap:
                if vMap[i] == cn:
                    vMap[i] = neighbor
            
            for i in vMap:
                if vMap[i] == neighbor:
                    counts[i] = tc
            
            cn_min_count = min(counts.values())
        
        for line in outlines:
            line[-1] = vMap[line[-1]]
        
        self.outlines = []
        for i, line in enumerate(outlines):
            if not len(self.outlines):
                self.outlines.append(line)
            else:
                if line[0] != self.outlines[-1][0]:
                    self.outlines.append(line)
                else:
                    if (line[3] == 0) and (i < len(outlines)-1) and (line[0] == outlines[i+1][0]) and (line[2] - line[1] < 5):
                        self.outlines[-1][2] = line[2]
                    else:
                        if line[3] != self.outlines[-1][3]:
                            self.outlines.append(line)
                        else:
                            self.outlines[-1][2] = line[2]
    
    def output(self, fil, pre='chr'):

        with open(fil, 'w') as out:
            for c, s, e, v in self.outlines:
                out.write('\t'.join([pre+c, str(s*self.binsize), str(e*self.binsize), str(v)])+'\n')

    
    def pieces(self, arr, minseq=10, scale='log'):

        idx = np.where(arr!=0)[0]
        tmp = np.split(idx, np.where(np.diff(idx)!=1)[0]+1)
        seqs = []
        for i in tmp:
            if i.size >= minseq:
                if scale == 'log':
                    seqs.append(np.log2(arr[i]))
                else:
                    seqs.append(arr[i])
        
        return seqs
    
    def assign_cnv(self, cnv_profile, hmmseq):

        seg = np.zeros(cnv_profile.size)
        for i in np.unique(hmmseq):
            idx = np.where(hmmseq == i)[0]
            tmp = np.split(idx, np.where(np.diff(idx)!=1)[0]+1)
            for t in tmp:
                local_mean = np.median(cnv_profile[t])
                seg[t] = local_mean
        
        return seg
    
    def call_intervals(self, cnvs):

        seg = []
        for i in np.unique(cnvs):
            idx = np.where(cnvs == i)[0]
            tmp = np.split(idx, np.where(np.diff(idx)!=1)[0]+1)
            for t in tmp:
                seg.append((t[0], t[-1]+1))
        
        seg.sort()

        return seg
    
    def get_states(self, maxdist=15):
        """
        Chromosome level estimation.
        """
        training_seqs = []
        whole_genome = np.r_[[]]
        for c in self.bin_cnv:
            arr = self.bin_cnv[c]
            nonzero = arr[arr > 0.1]
            idx = self.hampel_filter(np.log2(nonzero))
            filtered = nonzero[idx]
            whole_genome = np.r_[whole_genome, np.log2(filtered)]

            tmp = np.zeros(nonzero.size)
            tmp[idx] = filtered
            newarr = np.zeros(arr.size)
            newarr[arr > 0.1] = tmp

            self.bin_cnv[c] = newarr
            training_seqs.extend(self.pieces(newarr, scale='log'))
        
        X = whole_genome[:, np.newaxis]
        n_components_range = range(1, maxdist+1)
        _trace_gmm = {}
        aic_arr = []
        for n_components in n_components_range:
            gmm = GaussianMixture(n_components=n_components,
                                  covariance_type='diag',
                                  init_params='k-means++',
                                  n_init=10)
            gmm.fit(X)
            aic = gmm.aic(X)
            aic_arr.append(aic)
            _trace_gmm[n_components] = gmm
        
        aic_arr = np.r_[aic_arr]
        best_idx = aic_arr.argmin()
        gmm = _trace_gmm[best_idx+1]

        return training_seqs, gmm
    
    def hampel_filter(self, arr, ws = 10, n_sigmas = 3):
    
        n = len(arr)
        k = 1.4826 # scale factor for Gaussian distribution
        
        indices = []
        
        for i in range(n):
            if i < ws:
                si = 0
            else:
                si = i - ws
            if i + ws > n:
                ei = n
            else:
                ei = i + ws
            sub = arr[si:ei]
            x0 = np.median(sub)
            S0 = k * np.median(np.abs(sub - x0))
            
            if np.abs(arr[i] - x0) < n_sigmas * S0:
                indices.append(i)
        
        return indices
        