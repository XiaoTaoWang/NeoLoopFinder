import logging
import numpy as np
from neoloop.util import find_chrom_pre
from neoloop.callers import Fusion
from neoloop.util import map_coordinates, load_translocation_fil
from joblib import Parallel, delayed
import networkx as nx
from networkx.algorithms import all_shortest_paths, has_path

log = logging.getLogger(__name__)

def filterSV(clr, c1, c2, p1, p2, s1, s2, note, span, col, protocol, cutoff):

    fu = Fusion(clr, c1, c2, p1, p2, s1, s2, note, span=span, col=col, protocol=protocol, trim=False)
    fu.detect_bounds()
    fu.allele_slope()
    fu.correct_heterozygosity()
    if (fu.m_g_r > cutoff) and max(fu.u_s, fu.d_s) > 0.1:
        node = (c1, p1, s1, c2, p2, s2, note)
        Map = map_coordinates(fu.k_p, clr.binsize, c1, c2)[1]
        up = Map[fu.up_i]
        down = Map[fu.down_i]
        record = [node, (c1, up[1], fu.k_p[1]), (c2, fu.k_p[2], down[1])]
        return record
    else:
        return

def filterAssembly(clr, ID, span, col, path, protocol):

    wk = complexSV(clr, ID, span, col, protocol=protocol)
    wk.reorganize_matrix()
    wk.correct_heterozygosity()
    connectable = (len(wk.pair_binary)==wk.valid_counts) and wk.connectivity

    return connectable, path


class assembleSV(object):

    def __init__(self, clr, sv_fil, span=5000000, col='sweight', minIntra=500000,
        n_jobs=1, protocol='insitu', r_cutoff=0.6):

        self.clr = clr
        self.res = clr.binsize
        self.protocol = protocol
        self.span = span
        self.balance_type = col
        self.n_jobs = n_jobs
        # read SVs
        sv_list = load_translocation_fil(sv_fil, clr.binsize, minIntra)
        
        # filter sv by checking power-law decay of the induced interactions
        params = []
        for c1, p1, s1, c2, p2, s2, note in sv_list:
            params.append((clr, c1, c2, p1, p2, s1, s2, note, span, col, protocol, r_cutoff))
        
        results = Parallel(n_jobs=n_jobs, verbose=10)(delayed(filterSV)(*i) for i in params)
        self.queue = []
        self.lookup = {}
        for r in results:
            if not r is None:
                self.queue.append(r)
                self.lookup[r[0]] = [r[1], r[2]]
    
    def build_graph(self):

        nodes = {}
        for n, u, d in self.queue:
            k1 = (n, u, d)
            n2 = n[3:6] + n[:3] + (n[-1],)
            k2 = (n2, d, u)
            nodes[k1] = n
            nodes[k2] = n

        G = nx.DiGraph()
        for sv1 in nodes:
            for sv2 in nodes:
                if nodes[sv1] == nodes[sv2]:
                    continue
                dis = self.connect(sv1, sv2)
                if dis < np.inf:
                    G.add_weighted_edges_from([(sv1, sv2, dis)])
        
        self.graph = G
        self.nodemap = nodes
    
    def find_complexSV(self):
        
        pool = []
        nodes = self.graph.nodes()
        for n1 in nodes:
            for n2 in nodes:
                if n1 == n2:
                    continue
                pre = has_path(self.graph, n1, n2)
                if pre:
                    # do not consider edge weight
                    paths = list(all_shortest_paths(self.graph, n1, n2, weight=None))
                    for p in paths:
                        if not self._containloop(p):
                            pool.append((len(p), p))
        pool.sort(reverse=True)
        
        # so far, the candidate paths contain no self-loops, but are still redundant
        # check distance-decay for each pair of regions
        queue = [(self.clr, self._change_format(p[1]), self.span, self.balance_type, p[1], self.protocol) for p in pool]
        log.info('Filtering {0} redundant candidates ...'.format(len(queue)))
        jobs = Parallel(n_jobs=self.n_jobs, verbose=10)(delayed(filterAssembly)(*i) for i in queue)
        pre_alleles = []
        for ck, p in jobs:
            if ck:
                pre_alleles.append(p)
        
        # these assembly should exist within the same allele
        alleles = []
        for p in pre_alleles:
            for v in alleles:
                if self._issubset(p, v) or self._issubset(p, self._getreverse(v)):
                    break
            else:
                alleles.append(p)
        
        self.alleles = alleles
    
    def _containloop(self, path):

        loop_hit = []
        chain = []
        for p in path:
            chain.append(p[1])
            chain.append(p[2])
        
        for i in range(len(chain)):
            for j in range(len(chain)):
                if j - i > 1:
                    r1 = chain[i]
                    r2 = chain[j]
                    o, _ = self._check_overlap(r1, r2)
                    if o > 0.05:
                        loop_hit.append((r1, r2))
        
        return loop_hit
    
    def _getreverse(self, path):

        reverse = []
        for n, u, d in path[::-1]:
            n2 = n[3:6] + n[:3] + (n[-1],)
            k2 = (n2, d, u)
            reverse.append(k2)
        
        return reverse


    def _issubset(self, p1, p2):

        index = []
        for i in p1:
            if i in p2:
                index.append(p2.index(i))
            else:
                return False
        
        idx = np.r_[index]
        df = np.diff(idx)

        return np.all(df==1)
        
    
    def _check_overlap(self, i1, i2):

        if i1[0] != i2[0]:
            return 0, 0
        
        t1 = sorted(i1[1:])
        t2 = sorted(i2[1:])

        if (t1[1]<=t2[0]) or (t2[1]<=t1[0]):
            return 0, 0
        
        mi = t1 + t2
        mi.sort()
        OR = (mi[2] - mi[1]) / (mi[3] - mi[0]) # intersect / union

        return OR, mi[2]-mi[1]

    def connect(self, r1, r2):

        n1, u1, d1 = r1
        n2, u2, d2 = r2
        s11, s12 = n1[2], n1[5]
        s21, s22 = n2[2], n2[5]

        # check overlap between the 2nd part of sv1 and the 1st part of sv2
        ori = 0
        if s12 == '-':
            if s21 == '+':
                o, itr = self._check_overlap(d1, u2)
                if (o > ori) and (itr > 4*self.res):
                    ori = o
        else:
            if s21 == '-':
                o, itr = self._check_overlap(d1, u2)
                if (o > ori) and (itr > 4*self.res):
                    ori = o
        
        if ori > 0:
            return 1 / ori
        else:
            return np.inf
    
    def _change_format(self, path):

        names = []
        for p in path:
            tmp = ','.join(list(map(str, (p[0][-1],) + p[0][:-1])))
            names.append(tmp)
        
        b1 = path[0][1]
        b2 = path[-1][2]

        if path[0][0][2] == '+':
            t1 = [b1[0], str(min(b1[1:]))]
        else:
            t1 = [b1[0], str(max(b1[1:]))]
        
        if path[-1][0][5] == '-':
            t2 = [b2[0], str(max(b2[1:]))]
        else:
            t2 = [b2[0], str(min(b2[1:]))]
        
        t1 = ','.join(t1)
        t2 = ','.join(t2)

        line = names + [t1, t2]
        line = '\t'.join(line)
        
        return line
    
    def print(self, complexSVs, outfil):

        with open(outfil, 'w') as out:
            for path in complexSVs:
                line = self._change_format(path)
                out.write(line+'\n')
    
    def output_all(self, outfil):

        with open(outfil, 'w') as out:
            # alleles
            for i, path in enumerate(self.alleles):
                line = self._change_format(path)
                out.write('A{0}\t'.format(i)+line+'\n')
            # single SVs
            for i, sv in enumerate(self.lookup):
                path = [(sv, self.lookup[sv][0], self.lookup[sv][1])]
                line = self._change_format(path)
                out.write('C{0}\t'.format(i)+line+'\n')


class complexSV(Fusion):

    def __init__(self, clr, candidate, span=5000000, col='sweight', protocol='insitu',
        flexible=True, slopes={}):

        self.clr = clr
        self.res = clr.binsize
        self.protocol = protocol
        self.pre = find_chrom_pre(list(self.clr.chromnames))
        self.span = span
        self.flexible = flexible
        self.slopes = slopes
        
        if col in ['weight', 'sweight']:
            self.balance_type = col
        else:
            self.balance_type = False

        self.chroms = {}
        for c in self.clr.chromnames:
            self.chroms[c] = [0, self.clr.chromsizes[c]]
        
        self.parse_input(candidate)
        if len(self.sv_list)==1:
            self.tb, self.to = self.get_single_block(self.sv_list[0],
                                            left_bound=self.bounds[0],
                                            right_bound=self.bounds[1])
        else:
            self.tb = []
            self.to = []
            for i, sv in enumerate(self.sv_list):
                if i == 0:
                    intervals, directs = self.get_single_block(sv, left_bound=self.bounds[0])
                elif i == len(self.sv_list)-1:
                    intervals, directs = self.get_single_block(sv, right_bound=self.bounds[1])
                else:
                    intervals, directs = self.get_single_block(sv)
                self.tb.extend(intervals)
                self.to.extend(directs)
        
        self.load_expected()
    
    def parse_input(self, candidate):

        parse = candidate.split()
        self.sv_list = []
        for p in parse[:-2]:
            _, c1, p1, s1, c2, p2, s2 = p.split(',')
            p1, p2 = int(p1), int(p2)
            self.sv_list.append((c1, c2, p1, p2, s1, s2))

        b1 = parse[-2].split(',')
        b2 = parse[-1].split(',')
        self.bounds = [(b1[0], int(b1[1])), (b2[0], int(b2[1]))]

    
    def get_single_block(self, sv, left_bound=None, right_bound=None):

        c1, c2, p1, p2, s1, s2 = sv
        c1, c2 = self.pre+c1, self.pre+c2

        directions = [s1]
        if s2 == '+':
            directions.append('-')
        else:
            directions.append('+')
        
        intervals = []
        if s1 == '+':
            if self.flexible:
                if left_bound is None:
                    r1 = [c1, max(self.chroms[c1][0], p1 - self.span), p1]
                else:
                    r1 = [c1, left_bound[1], p1]
            else:
                r1 = [c1, max(self.chroms[c1][0], p1 - self.span), p1]
        else:
            if self.flexible:
                if left_bound is None:
                    r1 = [c1, p1, min(p1 + self.span, self.chroms[c1][1])]
                else:
                    r1 = [c1, p1, left_bound[1]]
            else:
                r1 = [c1, p1, min(p1 + self.span, self.chroms[c1][1])]
        intervals.append(r1)

        if s2 == '-':
            if self.flexible:
                if right_bound is None:
                    r2 = [c2, p2, min(p2 + self.span, self.chroms[c2][1])]
                else:
                    r2 = [c2, p2, right_bound[1]]
            else:
                r2 = [c2, p2, min(p2 + self.span, self.chroms[c2][1])]
        else:
            if self.flexible:
                if right_bound is None:
                    r2 = [c2, max(self.chroms[c2][0], p2 - self.span), p2]
                else:
                    r2 = [c2, right_bound[1], p2]
            else:
                r2 = [c2, max(self.chroms[c2][0], p2 - self.span), p2]
        intervals.append(r2)

        return intervals, directions
    
    def get_bias(self, r):

        if self.balance_type == False:
            col = 'weight'
        else:
            col = self.balance_type

        bias = self.clr.bins().fetch(r)[col].values
        bias[np.isnan(bias)] = 0

        return bias
    
    def get_matrix(self, r1, r2):

        M = self.clr.matrix(balance=self.balance_type).fetch(r1, r2)

        if self.balance_type:
            M[np.isnan(M)] = 0
        
        return M
    
    def reorganize_matrix(self):

        tb = self.tb
        to = self.to
        blocks = [tb[0]]
        orients = [to[0]]
        for i in range(1, len(tb)-1, 2):

            assert to[i]==to[i+1]
            assert tb[i][0]==tb[i+1][0]

            if to[i]=='+': # 5 --> 3
                b1 = tb[i][1]
                b2 = tb[i+1][2]
                assert b1 < b2
                blocks.append([tb[i][0], b1, b2])
                orients.append('+')
            else:
                b1 = tb[i][2]
                b2 = tb[i+1][1]
                assert b1 > b2
                blocks.append([tb[i][0], b2, b1])
                orients.append('-')
        blocks.append(tb[-1])
        orients.append(to[-1])

        # generate matrix and indexing
        index = []
        start = 0
        Map = []
        bounds = []
        for r, o in zip(blocks, orients):
            N = self.clr.bins().fetch(tuple(r))['chrom'].size
            index.append((start, start + N))
            coords = []
            for i in range(start, start + N):
                coords.append((r[0], r[1]//self.res*self.res + (i-start)*self.res))
            if o=='+':
                Map.extend(list(zip(coords, range(start, start + N))))
                bounds.append((start, (r[0], r[1]//self.res*self.res)))
                bounds.append((start + N - 1, (r[0], r[2]//self.res*self.res)))
            else:
                Map.extend(list(zip(coords[::-1], range(start, start + N))))
                bounds.append((start, (r[0], r[2]//self.res*self.res)))
                bounds.append((start + N - 1, (r[0], r[1]//self.res*self.res)))

            start = index[-1][1]

        forward = dict(Map)

        M = np.zeros((index[-1][1], index[-1][1]))
        bias_arr = np.r_[[]]
        for r1, o1, i1, j1 in zip(blocks, orients, index, range(len(blocks))):
            tmp = self.get_bias(r1)
            if o1=='+':
                bias_arr = np.r_[bias_arr, tmp]
            else:
                bias_arr = np.r_[bias_arr, tmp[::-1]]
            for r2, o2, i2, j2 in zip(blocks, orients, index, range(len(blocks))):
                if j1 > j2:
                    continue
                tmp = self.get_matrix(r1, r2)
                if (o1=='+') and (o2=='-'):
                    tmp = tmp[:,::-1]
                elif (o1=='-') and (o2=='+'):
                    tmp = tmp[::-1, :]
                elif (o1=='-') and (o2=='-'):
                    tmp = tmp[::-1, ::-1]
                M[i1[0]:i1[1], i2[0]:i2[1]] = tmp
        
        #x, y = M.nonzero()
        #M[y,x] = M[x,y]

        self.fusion_matrix = np.triu(M)
        self.Map = forward
        self.reverse = {forward[i]:i for i in forward}
        self.bounds = bounds
        self.orients = orients
        self.index = index
        self.blocks = blocks
        self.bias = bias_arr
        self.chains = {}
        for i, r in enumerate(self.index):
            for ri in range(r[0], r[1]):
                self.chains[ri] = i


    def correct_heterozygosity(self):

        Bias_M = self.bias[:,np.newaxis] * self.bias
        valid_bias = np.triu(Bias_M > 0)

        hcm = self.fusion_matrix.copy()
        slopes = {}
        exps = {}
        pair_binary = {}
        for i1, j1 in zip(self.index, range(len(self.blocks))):
            for i2, j2 in zip(self.index, range(len(self.blocks))):
                if j1 > j2:
                    continue
                inter_mask = np.zeros(hcm.shape, dtype=bool)
                inter_mask[i1[0]:i1[1], i2[0]:i2[1]] = True
                valid_pixel = inter_mask & valid_bias
                _rscores = []
                _slopes = []
                warning = []
                for maxdis in [400000, 1000000, 1500000, 2000000]:
                    local_exp = self._extract_xy(valid_pixel, mink=3, maxdis=maxdis)
                    r, s, w = self.linear_regression(local_exp, self.expected, min_samples=5)
                    warning.append(w)
                    _rscores.append(r)
                    _slopes.append(s)
                _rscores = np.r_[_rscores]
                r = _rscores.max()
                s = _slopes[_rscores.argmax()]
                    
                if any(warning):
                    pair_binary[(j1, j2)] = 0
                else:
                    pair_binary[(j1, j2)] = 1
                
                slopes[(j1, j2)] = [r, s]
                exps[(j1, j2)] = local_exp
        
        
        for i1, j1 in zip(self.index, range(len(self.blocks))):
            for i2, j2 in zip(self.index, range(len(self.blocks))):
                if j1 > j2:
                    continue
                ms = max(slopes[(j1, j1)][1], slopes[(j2, j2)][1])
                if (slopes[(j1, j2)][0] > 0.6) and (slopes[(j1, j2)][1] > 0) and (ms > 0):
                    if (j1, j2) in self.slopes:
                        hcm[i1[0]:i1[1], i2[0]:i2[1]] = hcm[i1[0]:i1[1], i2[0]:i2[1]] / self.slopes[(j1, j2)]
                    else:
                        if j1 == j2:
                            hcm[i1[0]:i1[1], i2[0]:i2[1]] = hcm[i1[0]:i1[1], i2[0]:i2[1]] / slopes[(j1, j2)][1]
                        else:
                            factor = min(1/slopes[(j1, j2)][1], 3/ms)
                            hcm[i1[0]:i1[1], i2[0]:i2[1]] = hcm[i1[0]:i1[1], i2[0]:i2[1]] * factor
        
        x, y = hcm.nonzero()
        hcm[y, x] = hcm[x, y]
        
        self.slopes_ = slopes
        self.exps = exps
        self.symm_hcm = hcm
        self.hcm = np.triu(hcm)
        valid_counts = sum([pair_binary[k] for k in pair_binary])
        self.valid_counts = valid_counts
        self.pair_binary = pair_binary
        # connectivity
        table = np.zeros((len(self.blocks), len(self.blocks)))
        for i, j in slopes:
            if slopes[i,j][0] > 0.6:
                table[i,j] = 1
                table[j,i] = 1
        connectivity = True
        for i in range(1, min(table.shape[0],3)):
            if np.any(table.diagonal(i)==0):
                connectivity = False
        self.connectivity = connectivity

    
    def index_to_coordinates(self, loop_list):

        loop_indices = []
        for i, j, v in loop_list:
            x = self.index_map[i]
            y = self.index_map[j]
            loop_indices.append((x, y, v))
        
        lines = {}
        for i, j, v in loop_indices:
            gdis = (j - i) * self.res
            loci1 = self.reverse[i]
            loci2 = self.reverse[j]
            if loci1 > loci2:
                loci1, loci2 = loci2, loci1
            key = (loci1[0], loci1[1], loci1[1]+self.res,
                   loci2[0], loci2[1], loci2[1]+self.res)
            if self.chains[i]==self.chains[j]:
                lines[key] = [gdis, 0] # normal loops
            else:
                lines[key] = [gdis, 1] # neo loops
        
        return lines