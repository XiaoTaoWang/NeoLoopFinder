#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import os, cooler, bisect, h5py
import numpy as np
from neoloop.util import find_chrom_pre
from cooler import util
from scipy.stats import rankdata, mode

# load CNV of any resolutions into a cool URI at specific resolution

class binCNV(object):

    def __init__(self, bedgraph, res):

        self.binsize = res

        ref_pre = ''

        self.load_bedGraph(bedgraph)
        ori_pre = find_chrom_pre(list(self.cnv_segment))
        D = {}
        for c in self.cnv_segment:
            ref_k = ref_pre+c.lstrip(ori_pre)
            chromlen = self.cnv_segment[c][-1][1]
            arr = []
            if self.cnv_res == res:
                for s, e, v in self.cnv_segment[c]:
                    arr.append(v)
            else:
                for i in range(0, chromlen, res):
                    start = i
                    end = min(i+res, chromlen)
                    tmp = self.calculate_bin(self.cnv_segment[c], [start, end])
                    arr.append(tmp)
            arr = np.r_[arr]
            D[ref_k] = arr
        
        # average copy numbers
        self.bin_cnv = D

    def calculate_bin(self, bychrom, interval):

        base_pairs = np.zeros(interval[1] - interval[0])
        idx = max(0, bisect.bisect(bychrom, interval)-1)
        for q in bychrom[idx:]:
            if q[1] <= interval[0]:
                continue
            if q[0] >= interval[1]:
                break
            s = q[0]-interval[0]
            if s < 0:
                s = 0
            e = q[1]-interval[0]
            base_pairs[s:e] += q[2]
        
        L = base_pairs.size
        if base_pairs[:L//2].mean() == base_pairs[L//2:].mean():
            return base_pairs.mean()
        else:
            return mode(base_pairs)[0][0]
    
    
    def load_bedGraph(self, infil):

        self.cnv_segment = {}
        with open(infil, 'r') as source:
            for line in source:
                p = line.rstrip().split()
                if p[0] in self.cnv_segment:
                    self.cnv_segment[p[0]].append([int(p[1]), int(p[2]), float(p[3])])
                else:
                    self.cnv_segment[p[0]] = [[int(p[1]), int(p[2]), float(p[3])]]
                    self.cnv_res = int(p[2]) - int(p[1])
        
        for c in self.cnv_segment:
            self.cnv_segment[c].sort()
        
    
    def assign_cnv(self, cooler_uri):
        
        cooler_lib = cooler.Cooler(cooler_uri)
        ref_pre = find_chrom_pre(cooler_lib.chromnames)
        cnv = np.r_[[]]
        for ref_k in cooler_lib.chromnames: # line with bin table
            bias = cooler_lib.bins().fetch(ref_k)['weight'].values
            c = ref_k.lstrip(ref_pre)
            if not c in self.bin_cnv:
                pre = np.zeros(len(bias))
                cnv = np.r_[cnv, pre]
                continue

            pre = self.bin_cnv[c]
            if len(bias) <= pre.size:
                pre = pre[:len(bias)]
            else:
                add = np.zeros(len(bias)-pre.size)
                pre = np.r_[pre, add]
            
            mask = np.isnan(bias) | (bias==0)
            pre[mask] = 0

            cnv = np.r_[cnv, pre]
        
        cnvi = rankdata(cnv, method='dense') - 1 # indices for quick Bias retrival
        
        # pre-check the CNV column
        cool_path, group_path = util.parse_cooler_uri(cooler_uri)
        with h5py.File(cool_path, 'r+') as h5:
            grp = h5[group_path]
            if 'CNV' in grp['bins']:
                del grp['bins']['CNV'] # Overwrite the CNV column
                del grp['bins']['CNVI']
        
        with h5py.File(cool_path, 'r+') as h5:
            grp = h5[group_path]
            h5opts = dict(compression='gzip', compression_opts=6)
            grp['bins'].create_dataset('CNV', data=cnv, **h5opts)
            grp['bins'].create_dataset('CNVI', data=cnvi, dtype='i1', **h5opts)

        del cnv, cnvi