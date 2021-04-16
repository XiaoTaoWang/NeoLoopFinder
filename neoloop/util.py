"""
Created on Tue Sep 18 22:15:12 2018

@author: XiaoTao Wang
"""
import numpy as np
import bisect, os

def find_chrom_pre(chromlabels):

    ini = chromlabels[0]
    if ini.startswith('chr'):
        return 'chr'
    
    else:
        return ''


def extract_sv_region(hic_pool, c1, c2, pos_1, pos_2, strand_1, strand_2,
                      radius=5000000):
    
    pre = find_chrom_pre(list(hic_pool.chromnames))
    c1 = pre+c1.lstrip('chr')
    c2 = pre+c2.lstrip('chr')

    res = hic_pool.binsize
    if c1==c2:
        if pos_1 > pos_2:
            pos_1, pos_2 = pos_2, pos_1
            strand_1, strand_2 = strand_2, strand_1
    
    p1 = pos_1 // res
    p2 = pos_2 // res
    r = radius // res

    chromsize_1 = hic_pool.chromsizes[c1]
    chromsize_2 = hic_pool.chromsizes[c2]

    if (strand_1=='+') and (strand_2=='+'):
        minx_1 = max(p1-r, 0)
        if c1!=c2:
            minx_2 = max(p2-r, 0)
        else:
            minx_2 = max(p2-r, p1+int((p2-p1)/2)) # handle short inversion ++
        k_p = [minx_1*res, min(p1*res+res,chromsize_1), min(p2*res+res,chromsize_2), minx_2*res]
    elif (strand_1=='+') and (strand_2=='-'):
        minx = max(p1-r, 0)
        maxx = min((p2+r+1)*res, chromsize_2)
        k_p = [minx*res, min(p1*res+res,chromsize_1), p2*res, maxx]
    elif (strand_1=='-') and (strand_2=='-'):
        if c1!=c2:
            maxx_1 = min((p1+r+1)*res, chromsize_1)
        else:
            maxx_1 = min((p1+r+1)*res, (p2-int((p2-p1)/2))*res) # handle short inversion --
        maxx_2 = min((p2+r+1)*res, chromsize_2)
        k_p = [maxx_1, p1*res, p2*res, maxx_2]
    else:
        maxx = min((p1+r+1)*res, chromsize_1)
        minx = max(p2-r, 0)
        k_p = [maxx, p1*res, min(p2*res+res,chromsize_2), minx*res]
    
    return k_p, pos_1, pos_2, strand_1, strand_2


def load_translocation_fil(fil_path, res, minIntra):

    lines = set()
    with open(fil_path, 'r') as source:
        for line in source:
            c1, c2, strands, i1, i2, note = line.rstrip().split()
            if '-' in i1:
                s, e = map(int, i1.split('-'))
                if e - s > res:
                    continue
                p1 = (s+e) // 2
            else:
                p1 = int(i1)
            
            if '-' in i2:
                s, e = map(int, i2.split('-'))
                if e - s > res:
                    continue
                p2 = (s+e) // 2
            else:
                p2 = int(i2)
            s1, s2 = strands[0], strands[1]
            # Filter by intra SV span
            if c1==c2:
                if abs(p2-p1) < minIntra:
                    continue
                if p1 > p2:
                    p1, p2 = p2, p1
                    s1, s2 = s2, s1
                if (s1 == '-') and (s2 == '+') and (abs(p2-p1) < minIntra*2): # duplicate
                    continue
                        
            lines.add((c1.lstrip('chr'), p1, s1, c2.lstrip('chr'), p2, s2, note))
    
    lines = sorted(lines)
    
    return lines

def parse_normal_peaks(peakfile):

    D = {}
    with open(peakfile, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            if len(parse)==3:
                chrom1 = chrom2 = parse[0].lstrip('chr')
                start1 = int(parse[1])
                start2 = int(parse[2])
            if len(parse)==6:
                chrom1, chrom2 = parse[0].lstrip('chr'), parse[3].lstrip('chr')
                start1 = int(parse[1])
                start2 = int(parse[4])
            if chrom1 > chrom2:
                chrom1, chrom2 = chrom2, chrom1
                start1, start2 = start2, start1
            if chrom1 == chrom2:
                if start1 > start2:
                    start1, start2 = start2, start1

            if (chrom1, chrom2) in D:
                D[(chrom1, chrom2)].append([start1, start2])
            else:
                D[(chrom1, chrom2)] = [[start1, start2]]
    for c in D:
        D[c].sort()
        D[c] = np.r_[D[c]]
    
    return D

def load_sv_peaks(peakfile, maxapart=10000000, onlyInter=True, last='res'):

    byregion = {}
    with open(peakfile, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            name = parse[0]
            res = int(parse[1])
            c1, p1, c2, p2 = parse[2:6]
            p1, p2 = int(p1), int(p2)
            apart = int(parse[6])
            if parse[-1] == 'False':
                induced = False
            else:
                induced = True
            if apart > maxapart:
                continue
            
            if onlyInter:
                if not induced:
                    continue

            if name in byregion:
                if last=='res':
                    byregion[name].append((c1, p1, c2, p2, int(res)))
                elif last=='full':
                    byregion[name].append((c1, p1, c2, p2, int(res), apart))
            else:
                if last=='res':
                    byregion[name] = [(c1, p1, c2, p2, int(res))]
                else:
                    byregion[name] = [(c1, p1, c2, p2, int(res), apart)]
    
    return byregion

def unique_peaks(byregion, max_mismatch=25000):

    from scipy.spatial import distance_matrix

    # byregion: returned by load_sv_peaks
    visited = {}
    peak_list = []
    for k in byregion:
        for peak in byregion[k]:
            c1, p1, c2, p2 = peak[:4]
            if c1 > c2:
                c1, c2 = c2, c1
                p1, p2 = p2, p1
            if not (c1, c2) in visited:
                visited[(c1, c2)] = np.r_[[[p1, p2]]]
                peak_list.append((c1, p1, c2, p2)+tuple(peak)[4:])
            else:
                dis = distance_matrix(visited[(c1,c2)], [[p1,p2]]).ravel()
                if dis.min() > max_mismatch:
                    visited[(c1, c2)] = np.r_[visited[(c1,c2)], [[p1,p2]]]
                    peak_list.append((c1, p1, c2, p2)+tuple(peak)[4:])
    
    return peak_list


def get_peak_loci(byregion, min_interval=15000):

    bychrom = {}
    for r in byregion:
        for i in byregion[r]:
            if i[4] >= min_interval:
                loci1 = [i[1], i[1]+i[4]]
                loci2 = [i[3], i[3]+i[4]]
            else:
                half = (min_interval - i[4]) // 2
                loci1 = [max(0,i[1]-half), i[1]+i[4]+half]
                loci2 = [max(0,i[3]-half), i[3]+i[4]+half]
            if not i[0] in bychrom:
                bychrom[i[0]] = [loci1]
            else:
                idx = max(0, bisect.bisect(bychrom[i[0]],loci1)-1)
                minp, maxp = loci1
                cache = []
                for q in bychrom[i[0]][idx:]:
                    if q[1] <= loci1[0]:
                        continue
                    if q[0] >= loci1[1]:
                        break
                    if q[0] < minp:
                        minp = q[0]
                    if q[1] > maxp:
                        maxp = q[1]
                    cache.append(q)
                for q in cache:
                    bychrom[i[0]].remove(q)
                bychrom[i[0]].append([minp, maxp])
            
            if not i[2] in bychrom:
                bychrom[i[2]] = [loci2]
            else:
                idx = max(0, bisect.bisect(bychrom[i[2]],loci2)-1)
                minp, maxp = loci2
                cache = []
                for q in bychrom[i[2]][idx:]:
                    if q[1] <= loci2[0]:
                        continue
                    if q[0] >= loci2[1]:
                        break
                    if q[0] < minp:
                        minp = q[0]
                    if q[1] > maxp:
                        maxp = q[1]
                    cache.append(q)
                for q in cache:
                    bychrom[i[2]].remove(q)
                bychrom[i[2]].append([minp, maxp])
            
            bychrom[i[0]].sort()
            bychrom[i[2]].sort()
    
    return bychrom

def set_match(peak_1, peak_2, max_mismatch=40000):

    from collections import defaultdict
    from scipy.spatial import distance_matrix

    only_1 = set()
    only_2 = set()
    overlap_1 = set()
    overlap_2 = set()
    visited = defaultdict(set)
    for r in peak_1:
        parse = r.split(',')
        c1, c2 = parse[1], parse[4]
        if not r in peak_2:
            for p in peak_1[r]:
                only_1.add((r,)+p)
            continue
        ref = [[i[1],i[3]] for i in peak_2[r]]
        query = [[i[1],i[3]] for i in peak_1[r]]
        for v in query:
            dis = distance_matrix([v],ref).ravel()
            argmin = dis.argmin()
            mindis = dis.min()
            if mindis<=min(0.2*abs(v[0]-v[1]), max_mismatch) and (not argmin in visited[r]): # tune the threshold if needed
                overlap_1.add((r,)+(c1,v[0],c2,v[1]))
                visited[r].add(argmin)
            else:
                only_1.add((r,)+(c1,v[0],c2,v[1]))
    
    for r in peak_2:
        for i in range(len(peak_2[r])):
            tmp = (r,) + peak_2[r][i]
            if i in visited[r]:
                overlap_2.add(tmp)
            else:
                only_2.add(tmp)
    
    return only_1, only_2, overlap_1, overlap_2


def peak_match(sv_peaks, normal_peaks, max_mismatch=100000):

    from collections import defaultdict
    from scipy.spatial import distance_matrix

    intra_sv_only = set()
    inter_sv_only = set()
    intra_normal_only = set()
    intra_sv_overlap = set()
    intra_normal_overlap = set()
    visited = defaultdict(set)
    for r in sv_peaks:
        parse = r.split(',')
        c1, c2 = parse[1], parse[4]
        if c1!=c2:
            for peak in sv_peaks[r]:
                inter_sv_only.add((r,)+peak)
            continue
        if not r in normal_peaks:
            for peak in sv_peaks[r]:
                intra_sv_only.add((r,)+peak)
            continue
        ref = [[i[1],i[3]] for i in normal_peaks[r]]
        query = [[i[1],i[3]] for i in sv_peaks[r]]
        for v in query:
            dis = distance_matrix([v],ref).ravel()
            argmin = dis.argmin()
            mindis = dis.min()
            if mindis<=min(0.2*abs(v[0]-v[1]), max_mismatch) and (not argmin in visited[r]): # tune the threshold if needed
                intra_sv_overlap.add((r,)+(c1,v[0],c2,v[1]))
                visited[r].add(argmin)
            else:
                intra_sv_only.add((r,)+(c1,v[0],c2,v[1]))
    
    for r in normal_peaks:
        for i in range(len(normal_peaks[r])):
            tmp = (r,) + normal_peaks[r][i]
            if i in visited[r]:
                intra_normal_overlap.add(tmp)
            else:
                intra_normal_only.add(tmp)          
    
    return intra_sv_only, intra_sv_overlap, intra_normal_only, intra_normal_overlap, inter_sv_only

def map_coordinates(k_p, res, chrom_1, chrom_2, unit='bp'):
    
    # Coordinates on chrom 1
    if k_p[0] < k_p[1]:
        endbin = k_p[1]//res
        if endbin*res < k_p[1]:
            endbin += 1
        startbin = k_p[0]//res
        pl_1 = endbin - startbin
        if unit=='bp':
            D1 = dict(zip(range(k_p[0],k_p[1],res), range(0,pl_1)))
        else:
            D1 = dict(zip(range(startbin,endbin), range(0,pl_1)))
    else:
        endbin = k_p[0]//res
        if endbin*res < k_p[0]:
            endbin += 1
        startbin = k_p[1]//res
        pl_1 = endbin - startbin
        if unit=='bp':
            D1 = dict(zip(range(k_p[1],k_p[0],res), range(pl_1-1,-1,-1)))
        else:
            D1 = dict(zip(range(startbin,endbin), range(pl_1-1,-1,-1)))
    
    # Coordinates on chrom 2
    if k_p[2] < k_p[3]:
        endbin = k_p[3]//res
        if endbin*res < k_p[3]:
            endbin += 1
        startbin = k_p[2]//res
        pl_2 = endbin - startbin
        if unit=='bp':
            D2 = dict(zip(range(k_p[2],k_p[3],res), range(pl_1,pl_1+pl_2)))
        else:
            D2 = dict(zip(range(startbin,endbin), range(pl_1,pl_1+pl_2)))
    else:
        endbin = k_p[2]//res
        if endbin*res < k_p[2]:
            endbin += 1
        startbin = k_p[3]//res
        pl_2 = endbin - startbin
        if unit=='bp':
            D2 = dict(zip(range(k_p[3],k_p[2],res), range(pl_1+pl_2-1,pl_1-1,-1)))
        else:
            D2 = dict(zip(range(startbin,endbin), range(pl_1+pl_2-1,pl_1-1,-1)))
    
    forward = {} # global coordinates to region bin index
    for k in D1:
        forward[(chrom_1, k)] = D1[k]
    for k in D2:
        forward[(chrom_2, k)] = D2[k]
    
    reverse = {v:k for k,v in forward.items()}
    
    return forward, reverse

def properU(pos):
    
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])
