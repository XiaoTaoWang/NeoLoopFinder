#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import logging, copy
from collections import defaultdict
import numpy as np
from sklearn.mixture import GaussianMixture
from neoloop.cnv.loadcnv import binCNV
from pomegranate import NormalDistribution, HiddenMarkovModel, GeneralMixtureModel, State

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
    Recursively segment the interval x[start:end] returning a list L of pairs (i,j) where each (i,j) is a significant segment.
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

def combine_segment(sig, cbs_seg, hmm_seg, bufsize=4):

    start = cbs_seg[0][0]
    pool1 = set()
    for b in cbs_seg:
        for i in range(bufsize):
            pool1.add(b[1]-i)
            pool1.add(b[1]+i)

    pool2 = set([b[1] for b in hmm_seg])
    pool = [start] + sorted(pool1 & pool2)
    seg = []
    for i in range(len(pool)-1):
        s, e = pool[i], pool[i+1]
        v = np.median(sig[s:e])
        seg.extend([v]*(e-s))
    seg = np.r_[seg]

    return seg

def impute(self, arr):

    zeros = np.where(arr==0)[0]
    nonzeros = np.where(arr!=0)[0]
    fill = arr.copy()
    for i in zeros:
        idx = nonzeros.searchsorted(i)
        if (idx == 0) or (idx == nonzeros.size):
            continue # leading or tailing 0s
        ri = nonzeros[idx]
        li = nonzeros[idx-1]
        if arr[li] < arr[ri]:
            fill[i] = arr[li]
        else:
            fill[i] = arr[ri]
    
    return fill

class HMMsegment(binCNV):

    def __init__(self, bedgraph, res, nproc=1, ploidy=2):

        binCNV.__init__(self, bedgraph, res)
        self.n_jobs = nproc
        self.ploidy = ploidy

    def segment(self, min_seg=20):

        outlines = []
        counts = defaultdict(int)
        vMap = {}
        self._original = {}
        for c in self.bin_cnv:
            logger.info('Segmenting Chromosome {0} ...'.format(c))
            sig = self.bin_cnv[c]
            sig, hmm_seg, scale = self._segment(sig) # pure HMM-based segmentation
            self._original[c] = [sig, hmm_seg]
            logger.info('Combine results from the CBS algorithm ...')
            if scale=='log':
                logsig = np.zeros_like(sig)
                logsig[sig>0] = np.log2(sig[sig>0])
                cbs_seg = segment(logsig)
            else:
                cbs_seg = segment(sig)
            seg = combine_segment(sig, cbs_seg, hmm_seg) # combine results from cbs and hmm
            
            cur = [c, 0, 1, seg[0]]
            for i in range(1, seg.size):
                v = seg[i]
                if v == cur[-1]:
                    cur[2] = cur[2] + 1
                else:
                    tmp = sig[cur[1]:cur[2]]
                    if tmp.sum() == 0:
                        cn = 0
                    else:
                        cn = int(np.round(np.median(tmp[tmp!=0])*self.ploidy))
                    cur[3] = cn
                    vMap[cn] = cn
                    counts[cn] += (cur[2]-cur[1])
                    outlines.append(cur)
                    start = cur[2]
                    cur = [c, start, start+1, v]
            
            tmp = sig[cur[1]:cur[2]]
            if tmp.sum() == 0:
                cn = 0
            else:
                cn = int(np.round(np.median(tmp[tmp!=0])*self.ploidy))
            cur[3] = cn
            vMap[cn] = cn
            counts[cn] += (cur[2]-cur[1])
            outlines.append(cur)
        
        cn_min_count = min(counts.values())
        while cn_min_count < min_seg:
            for i in counts:
                if counts[i]==cn_min_count:
                    cn = i
            neighbor = 0
            mindis = 1000 # a very large number to start
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

    def _segment(self, arr, components=2):

        nonzero = arr[arr > 0]
        idx = self.hampel_filter(np.log2(nonzero))
        filtered = nonzero[idx]

        log_gmm = self.get_states(np.log2(filtered))
        log_means, log_probs = log_gmm.means_.ravel(), log_gmm.weights_
        ln_gmm = self.get_states(filtered) # to improve the sensitivity
        ln_means, ln_probs = ln_gmm.means_.ravel(), ln_gmm.weights_
        if (len(log_means) == 1):
            means, probs = ln_means, ln_probs
            scale = 'linear'
        else:
            means, probs = log_means, log_probs
            scale = 'log'

        logger.info('Estimated HMM state number: {0} ({1} scale)'.format(len(means), scale))
        model = HiddenMarkovModel()
        # GMM emissions
        dists = []
        for m in means:
            tmp = []
            for i in range(components):
                e = m + (-1)**i * ((i+1)//2) * 0.5
                s = 0.5
                tmp.append(NormalDistribution(e, s))
            mixture = State(GeneralMixtureModel(tmp), name=str(m))
            dists.append(mixture)
        model.add_states(*tuple(dists))
        # transition matrix
        for i in range(len(means)):
            for j in range(len(means)):
                if i==j:
                    model.add_transition(dists[i], dists[j], 0.8)
                else:
                    model.add_transition(dists[i], dists[j], 0.2/(len(means)-1))
        
        # starts and ends
        for i in range(len(means)):
            model.add_transition(model.start, dists[i], probs[i])
        
        model.bake()

        # training sequences
        tmp = np.zeros(nonzero.size)
        tmp[idx] = filtered
        newarr = np.zeros(arr.size)
        newarr[arr > 0] = tmp

        if len(means) > 1:
            model.fit(self.pieces(newarr, scale=scale), algorithm='baum-welch', n_jobs=self.n_jobs,
                    max_iterations=5000, stop_threshold=2e-4)
            
            queue = newarr[newarr > 0]
            
            if scale=='log':
                seq = np.r_[[s.name for i, s in model.viterbi(np.log2(queue))[1][1:]]]
            else:
                seq = np.r_[[s.name for i, s in model.viterbi(queue)[1][1:]]]
            seg = self.assign_cnv(queue, seq)
            
            predicted = np.zeros(newarr.size)
            predicted[newarr > 0] = seg
            seg = self.call_intervals(predicted)
        else:
            seg = [(0, newarr.size)]
        
        return newarr, seg, scale
    
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
    
    def get_states(self, arr, maxdist=15):
        """
        Chromosome level estimation.
        """
        X = arr[:, np.newaxis]

        n_components_range = range(1, maxdist+1)
        _trace_gmm = {}
        bic_arr = []
        for n_components in n_components_range:
            gmm = GaussianMixture(n_components=n_components,
                                covariance_type='diag')
            gmm.fit(X)
            bic = gmm.bic(X)
            bic_arr.append(bic)
            _trace_gmm[n_components] = gmm
        
        bic_arr = np.r_[bic_arr]
        best_idx = bic_arr.argmin()
        #loss_rate = np.abs((bic_arr - bic_arr.min()) / bic_arr.min())
        #best_idx = np.where(loss_rate < 0.01)[0][-1]
        
        gmm = _trace_gmm[best_idx+1]
        
        return gmm

    
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
        