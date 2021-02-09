#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

from neoloop.callers import Fusion
from neoloop.cnv.loadcnv import binCNV
from neoloop.util import find_chrom_pre
import numpy as np
from joblib import Parallel, delayed
import os

def outside_trans(arg, **kwarg):

    return Train.trans_sum(*arg, **kwarg)


def outside_cis(arg, **kwarg):

    return Train.cis_vec(*arg, **kwarg)


class Train(binCNV, Fusion):

    def __init__(self, clr, cnv_file, n_jobs = 1):
        
        self.clr = clr
        self.res = self.clr.binsize
        self.n_jobs = n_jobs

        binCNV.__init__(self, cnv_file, self.res)

        self.queue = {}
        self.pre = find_chrom_pre(clr.chromnames)
        for c in self.bin_cnv:
            tmp = self.seg(c)
            for v in tmp:
                if not v in self.queue:
                    self.queue[v] = []
                self.queue[v].extend(tmp[v])
    
    def expected_by_cnv(self):

        # Parallel trans stats computation
        trans = {}
        for ci in self.queue:
            for cj in self.queue:
                if ci > cj:
                    continue
                key = (ci, cj)
                params = set()
                for r1 in self.queue[ci]:
                    for r2 in self.queue[cj]:
                        rr1, rr2 = r1, r2
                        if rr1 > rr2:
                            rr1, rr2 = rr2, rr1
                        if rr1[0] == rr2[0]:
                            continue
                        params.add((self, rr1, rr2))
                if len(params):
                    results = Parallel(n_jobs=self.n_jobs)(delayed(outside_trans)(i) for i in params)
                    results = np.r_[results]
                    counts = results[:,0].sum()
                    num = results[:,1].sum()
                    if num > 0:
                        trans[key] = [counts / num, num]
                    else:
                        trans[key] = [0, 0]
                else:
                    trans[key] = [0, 0]

        self.trans = trans
        
        # Parallel cis stats computation
        cis = {}
        for ci in self.queue:
            for cj in self.queue:
                if ci > cj:
                    continue
                key = (ci, cj)
                params = set()
                for r1 in self.queue[ci]:
                    for r2 in self.queue[cj]:
                        rr1, rr2 = r1, r2
                        if rr1 > rr2:
                            rr1, rr2 = rr2, rr1
                        if rr1[0] != rr2[0]:
                            continue
                        params.add((self, rr1, rr2))
                tmp = {}
                if len(params):
                    results = Parallel(n_jobs=self.n_jobs)(delayed(outside_cis)(i) for i in params)
                    for D in results:
                        for i in D:
                            if not i in tmp:
                                tmp[i] = [0, 0]
                            tmp[i][0] += D[i][0]
                            tmp[i][1] += D[i][1]
                cis[key] = {}
                for i in tmp:
                    if tmp[i][1] >= 3:
                        cis[key][i] = [tmp[i][0] / tmp[i][1], tmp[i][1]]
        
        self.cis = cis
    
    def trans_factor(self, ref = 3):

        tmp = {}
        ref_v = self.trans[(ref, ref)][0]
        for i in self.trans:
            f1 = self.trans[i][0] / ref_v
            tmp[i] = f1
        
        return tmp
    
    def cis_factor(self, ref = 3):

        slopes = {}
        intervals = [(0, 20000000//self.res),
                     (20000000//self.res, 500000000//self.res)]
        
        k = (ref, ref)
        ref_e = self.cis[k]
        for i in sorted(self.cis):
            slopes[i] = [[1, 1]]
            if i == k:
                slopes[i] = [[1, 1]]
            else:
                cur = self.cis[i]
                pool = []
                for minx, maxx in intervals:
                    exp2 = {t:ref_e[t][0] for t in ref_e if (t >= minx) and (t <= maxx)}
                    exp1 = {t:cur[t][0] for t in cur if (t >= minx) and (t <= maxx)}
                    cr, cs, w = self.linear_regression(exp1, exp2)
                    pool.append([cs, cr])
                slopes[i] = pool

        return slopes
    
    def seg(self, chrom, minsize=2):

        arr = self.bin_cnv[chrom]
        chrom = self.pre + chrom
        chromLen = self.clr.chromsizes[chrom]
        values = np.unique(arr)
        D = {} # cnv-centered segments
        for v in values:
            tmp = np.where(arr==v)[0]
            pieces = np.split(tmp, np.where(np.diff(tmp)!=1)[0]+1)
            filtered = [p for p in pieces if len(p) >= minsize]
            intervals = []
            for p in filtered:
                intervals.append((chrom, p[0]*self.res, min((p[-1]+1)*self.res, chromLen)))
            
            if len(intervals):
                key = float(v)
                D[key] = intervals
        
        return D
    
    def trans_sum(self, r1, r2):

        M = self.clr.matrix(balance=False, sparse=True).fetch(r1, r2).tocsr()
        M2 = self.clr.matrix(balance=True, sparse=False).fetch(r1, r2)
        x, y = np.where(np.isfinite(M2))

        if x.size > 0:
            total = M[x, y].sum()
            num = x.size
            return [total, num]
        else:
            return [0, 0]
        
    
    def cis_vec(self, r1, r2):
        '''
        r1 and r2 must be on the same chromosome, and r2 > r1
        '''
        M = self.clr.matrix(balance=False, sparse=True).fetch(r1, r2).tocsr()
        M2 = self.clr.matrix(balance=True, sparse=False).fetch(r1, r2)
        x, y = np.where(np.isfinite(M2))
        D = {}
        if x.size > 0:
            gx = x + r1[1]//self.res
            gy = y + r2[1]//self.res

            dis = gy - gx
            maxi = dis.max()
            mini = max(1, dis.min())
            for i in range(mini, maxi+1):
                mask = dis == i
                xi, yi = x[mask], y[mask]
                if xi.size > 0:
                    D[i] = [M[xi, yi].sum(), xi.size]
        
        return D


class SimCNV(Train):

    def __init__(self, clr, cnv_file, trans_factor={}, cis_factor={}, out_folder='SimCNV'):
        
        Train.__init__(self, clr, cnv_file)

        self.trans = trans_factor
        self.cis = cis_factor
        self.out_folder = out_folder
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)
        self.cnv_by_chrom()
    
    def cnv_by_chrom(self):

        pairs = {}
        for i, c1 in enumerate(self.clr.chromnames):
            for j, c2 in enumerate(self.clr.chromnames):
                if i > j:
                    continue
                tmp = set()
                for k1 in self.queue:
                    for k2 in self.queue:
                        if k1 > k2:
                            continue
                        key = (k1, k2)
                        for r1 in self.queue[k1]:
                            for r2 in self.queue[k2]:
                                if (r1[0]==c1) and (r2[0]==c2):
                                    tmp.add((r1, r2, key))
                                elif (r2[0]==c1) and (r1[0]==c2):
                                    tmp.add((r2, r1, key))

                pairs[(c1, c2)] = tmp
        
        self.pairs = pairs
                        

    def correct_trans(self, chrom1, chrom2, ref=3):

        M = self.clr.matrix(balance=False).fetch(chrom1, chrom2)
        pairs = self.pairs[(chrom1, chrom2)]
        for r1, r2, k in pairs:
            if k==(ref, ref):
                continue
            s1, e1 = r1[1]//self.res, r1[2]//self.res
            s2, e2 = r2[1]//self.res, r2[2]//self.res
            f = self.trans[k]
            M[s1:e1, s2:e2] = M[s1:e1, s2:e2] * f
        
        self.out_matrix(M, chrom1, chrom2)
    
    def correct_cis(self, chrom, ref=3):

        M = self.clr.matrix(balance=False).fetch(chrom)
        pairs = self.pairs[(chrom, chrom)]

        intervals = [(0, 20000000//self.res),
                     (20000000//self.res, 500000000//self.res)]

        for r1, r2, k in pairs:
            if r1 > r2:
                continue
            s1, e1 = r1[1]//self.res, r1[2]//self.res
            s2, e2 = r2[1]//self.res, r2[2]//self.res
            fs = self.cis[k]
            if len(fs)!=len(intervals):
                continue

            f1, rsquare1 = fs[0]
            f2, rsquare2 = fs[1]
            if f1 > 0 and rsquare1 > 0.5:
                f = f1
            elif f2 > 0 and rsquare2 > 0.5:
                f = f2
            else:
                f = 1
            
            M[s1:e1, s2:e2] = M[s1:e1, s2:e2] * f
                    

        self.out_matrix(M, chrom, chrom)


    def out_matrix(self, M, chrom1, chrom2):

        outfil = os.path.join(self.out_folder, '{0}_{1}.txt'.format(chrom1.lstrip('chr'), chrom2.lstrip('chr')))
        if chrom1 == chrom2:
            M = np.triu(M)
        x, y = M.nonzero()
        values = M[x, y]
        with open(outfil, 'w') as out:
            for i, j, v in zip(x, y, values):
                line = '{0}\t{1}\t{2}\n'.format(i, j, v)
                out.write(line)
            
    