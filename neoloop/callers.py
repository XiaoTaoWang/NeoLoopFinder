"""
Created on Wed Jun 27 18:48:42 2018

@author: XiaoTao Wang
"""

import logging, neoloop, os, joblib, math
import numpy as np
from scipy import sparse
from scipy.stats import poisson
from sklearn.linear_model import HuberRegressor
from sklearn.isotonic import IsotonicRegression
from neoloop.util import find_chrom_pre, map_coordinates, distance_normalize, image_normalize
from scipy.stats import spearmanr
from scipy.ndimage import gaussian_filter

log = logging.getLogger(__name__)

def check_increasing(x, y):
    """Determine whether y is monotonically correlated with x.
    y is found increasing or decreasing with respect to x based on a Spearman
    correlation test.
    Parameters
    ----------
    x : array-like of shape (n_samples,)
            Training data.
    y : array-like of shape (n_samples,)
        Training target.
    Returns
    -------
    increasing_bool : boolean
        Whether the relationship is increasing or decreasing.
    Notes
    -----
    The Spearman correlation coefficient is estimated from the data, and the
    sign of the resulting estimate is used as the result.
    In the event that the 95% confidence interval based on Fisher transform
    spans zero, a warning is raised.
    References
    ----------
    Fisher transformation. Wikipedia.
    https://en.wikipedia.org/wiki/Fisher_transformation
    """

    # Calculate Spearman rho estimate and set return accordingly.
    rho, _ = spearmanr(x, y)
    warning = False
    increasing_bool = rho >= 0

    # Run Fisher transform to get the rho CI, but handle rho=+/-1
    if rho not in [-1.0, 1.0] and len(x) > 3:
        F = 0.5 * math.log((1. + rho) / (1. - rho))
        F_se = 1 / math.sqrt(len(x) - 3)

        # Use a 95% CI, i.e., +/-1.96 S.E.
        # https://en.wikipedia.org/wiki/Fisher_transformation
        rho_0 = math.tanh(F - 1.96 * F_se)
        rho_1 = math.tanh(F + 1.96 * F_se)

        # Warn if the CI spans zero.
        if np.sign(rho_0) != np.sign(rho_1):
            warning = True

    return warning, increasing_bool

def clean_inputs(Xi, X, Y, increasing_bool):

    IR = IsotonicRegression(increasing=increasing_bool)
    IR.fit(Xi, Y)
    Y1 = IR.predict(Xi)
    vi = np.where(np.diff(Y1) < 0)[0]
    pieces = np.split(vi, np.where(np.diff(vi)!=1)[0]+1)
    si = 0
    for i in range(len(pieces)-1):
        p1 = pieces[i]
        p2 = pieces[i+1]
        if p1.size / (p2[0] - p1[0]) > 0.5:
            si = p1[0]
            break
    
    if si / X.size > 0.3: # if more than 1/4 data discarded
        si = vi[0]
        if si / X.size > 0.3:
            si = 0
    
    X = X[si:]
    Y = Y[si:]

    return X, Y


class Fusion(object):

    def __init__(self, clr, c1, c2, p1, p2, s1, s2, sv_type, span=5000000,
        col='sweight', trim=True, protocol='insitu', expected_values=None):

        self.clr = clr
        self.p1 = p1
        self.p2 = p2
        self.s1 = s1
        self.s2 = s2
        self.res = clr.binsize
        pre = find_chrom_pre(list(self.clr.chromnames))
        self.c1 = pre + c1.lstrip('chr')
        self.c2 = pre + c2.lstrip('chr')
        self.chromsize1 = self.clr.chromsizes[self.c1]
        self.chromsize2 = self.clr.chromsizes[self.c2]
        self.balance_type = col
        self.protocol = protocol

        self.name = ','.join(map(str,
                                [sv_type, c1.lstrip('chr'), p1, s1, c2.lstrip('chr'), p2, s2]))

        self.get_matrix(span, col, trim)
        if expected_values is None:
            self.load_expected()
        else:
            self.expected = expected_values
    
    def load_expected(self):

        data_folder = os.path.join(os.path.split(neoloop.__file__)[0], 'data')
        paths = {
            5000: os.path.join(data_folder, 'gm.insitu.expected.5k.pkl'),
            10000: os.path.join(data_folder, 'gm.insitu.expected.10k.pkl'),
            20000: os.path.join(data_folder, 'gm.insitu.expected.20k.pkl'),
            25000: os.path.join(data_folder, 'gm.insitu.expected.25k.pkl'),
            40000: os.path.join(data_folder, 'gm.insitu.expected.40k.pkl'),
            50000: os.path.join(data_folder, 'gm.insitu.expected.50k.pkl')
        }
        
        self.expected = {}
        if self.res in paths:
            self.expected = joblib.load(paths[self.res])
         

    def get_matrix(self, span, col, trim):

        res = self.res
        chromsize1 = self.chromsize1
        chromsize2 = self.chromsize2
        p1 = self.p1 // res
        p2 = self.p2 // res
        r = span // res
        
        # 4 kind of directions
        if (self.s1 == '+') and (self.s2 == '+'):
            minx_1 = max(p1-r, 0)
            if self.c1 != self.c2:
                minx_2 = max(p2 - r, 0)
            else:
                if trim:
                    minx_2 = max(p2 - r, p1 + int((p2 - p1) / 2)) # handle short inversion ++
                else:
                    minx_2 = max(p2 - r, p1 + 1)

            k_p = [minx_1 * res, min(p1 * res + res, chromsize1),
                   min(p2 * res + res, chromsize2), minx_2 * res]
            # part 1
            M1 = self.clr.matrix(balance=col).fetch((self.c1, k_p[0], k_p[1]))
            # part 2
            M2 = self.clr.matrix(balance=col).fetch((self.c2, k_p[3], k_p[2]))
            M2 = M2[::-1,::-1]
            # part 3
            M3 = self.clr.matrix(balance=col).fetch((self.c1, k_p[0], k_p[1]),
                                                    (self.c2, k_p[3], k_p[2]))
            M3 = M3[:,::-1]
        elif (self.s1 == '+') and (self.s2 == '-'):
            minx = max(p1 - r, 0)
            maxx = min((p2 + r + 1) * res, chromsize2)
            k_p = [minx * res, min(p1 * res + res, chromsize1),
                   p2 * res, maxx]
            # part 1
            M1 = self.clr.matrix(balance=col).fetch((self.c1, k_p[0], k_p[1]))
            # part 2
            M2 = self.clr.matrix(balance=col).fetch((self.c2, k_p[2], k_p[3]))
            # part 3
            M3 = self.clr.matrix(balance=col).fetch((self.c1, k_p[0], k_p[1]),
                                                    (self.c2, k_p[2], k_p[3]))
        elif (self.s1 == '-') and (self.s2 == '-'):
            if self.c1 != self.c2:
                maxx_1 = min((p1 + r + 1) * res, chromsize1)
            else:
                if trim:
                    maxx_1 = min((p1 + r + 1) * res,
                                 (p2 - int((p2 - p1) / 2)) * res) # handle short inversion --
                else:
                    maxx_1 = min((p1 + r + 1) * res, p2 * res - res)
            maxx_2 = min((p2 + r + 1) * res, chromsize2)
            k_p = [maxx_1, p1 * res, p2 * res, maxx_2]
            # part 1
            M1 = self.clr.matrix(balance=col).fetch((self.c1, k_p[1], k_p[0]))
            M1 = M1[::-1,::-1] # - --> +
            # part 2
            M2 = self.clr.matrix(balance=col).fetch((self.c2, k_p[2], k_p[3]))
            # part 3
            M3 = self.clr.matrix(balance=col).fetch((self.c1, k_p[1], k_p[0]),
                                                    (self.c2, k_p[2], k_p[3]))
            M3 = M3[::-1,:] # - --> +, + --> -
        else:
            if self.c1 != self.c2:
                maxx = min((p1 + r + 1) * res, chromsize1)
                minx = max(p2 - r, 0)
                k_p = [maxx, p1 * res,
                    min(p2 * res + res, chromsize2), minx * res] # no intra SV is -+ ?
            else:
                maxx_1 = min(p1 + r + 1, p2 - int((p2 - p1) / 2) - 1)
                minx_2 = max(p2 - r, p1 + int((p2 - p1) / 2) + 1)
                k_p = [maxx_1 * res, p1 * res, min(p2 * res + res, chromsize2), minx_2 * res]

            # part 1
            M1 = self.clr.matrix(balance=col).fetch((self.c1, k_p[1], k_p[0]))
            M1 = M1[::-1,::-1] # - --> +
            # part 2
            M2 = self.clr.matrix(balance=col).fetch((self.c2, k_p[3], k_p[2]))
            M2 = M2[::-1,::-1] # - --> +
            # part 3
            M3 = self.clr.matrix(balance=col).fetch((self.c1, k_p[1], k_p[0]),
                                                    (self.c2, k_p[3], k_p[2]))
            M3 = M3[::-1,::-1] # - --> +, - --> +
        
        # put three parts together
        size = M1.shape[0] + M2.shape[0]
        M = np.zeros((size, size))
        M[:M1.shape[0], :M1.shape[0]] = M1
        M[M1.shape[0]:, M1.shape[0]:] = M2
        M[:M1.shape[0], M1.shape[0]:] = M3
        M = np.triu(M)
    
        p1L = M1.shape[0]
        
        if col:
            M[np.isnan(M)] = 0
        
        self.fusion_matrix = M
        self.k_p = k_p
        self.fusion_point = p1L

    
    def extract_bias(self, col='sweight'):
        
        k_p = self.k_p
        if k_p[0] < k_p[1]:
            t1 = self.clr.bins().fetch((self.c1, k_p[0], k_p[1]))[col].values
        else:
            t1 = self.clr.bins().fetch((self.c1, k_p[1], k_p[0]))[col].values
            t1 = t1[::-1]
        
        if k_p[2] < k_p[3]:
            t2 = self.clr.bins().fetch((self.c2, k_p[2], k_p[3]))[col].values
        else:
            t2 = self.clr.bins().fetch((self.c2, k_p[3], k_p[2]))[col].values
            t2 = t2[::-1]
        
        m1 = np.logical_not((t1==0) | np.isnan(t1))
        w1 = np.zeros_like(t1)
        w1[m1] = 1 / t1[m1]
        m2 = np.logical_not((t2==0) | np.isnan(t2))
        w2 = np.zeros_like(t2)
        w2[m2] = 1 / t2[m2]
    
        self.bias = np.r_[w1, w2]


    def detect_bounds(self):

        from sklearn.decomposition import PCA
        from scipy.ndimage import gaussian_filter

        n = self.fusion_matrix.shape[0]
        junc = self.fusion_point

        # Detect bounds by using the induced contacts
        new = self.fusion_matrix.copy()
        up_i = 0
        down_i = n - 1
        up_var_ratio = 0
        down_var_ratio = 0
        # locate the upstream and downstream bound independently
        inter = self.fusion_matrix[:junc, junc:]
        rowmask = inter.sum(axis=1) != 0
        colmask = inter.sum(axis=0) != 0
        if (rowmask.sum() >= 10) and (colmask.sum() >= 10):
            # row, upstream
            new = inter[rowmask][:,colmask]
            corr = gaussian_filter(np.corrcoef(new, rowvar=True), sigma=1)
            try:
                pca = PCA(n_components=3, whiten=True)
                pca1 = pca.fit_transform(corr)[:,0]
                up_var_ratio = pca.explained_variance_ratio_[0]
                loc = self.locate(pca1, left_most=False)
                up_i = np.where(rowmask)[0][loc] # included
            except:
                up_i = 0
                up_var_ratio = 0

            # column, downstream
            corr = gaussian_filter(np.corrcoef(new, rowvar=False), sigma=1)
            try:
                pca = PCA(n_components=3, whiten=True)
                pca1 = pca.fit_transform(corr)[:,0]
                down_var_ratio = pca.explained_variance_ratio_[0]
                loc = self.locate(pca1, left_most=True)
                down_i = np.where(colmask)[0][loc] + junc
            except:
                down_i = n - 1
                down_var_ratio = 0
        
        self.up_i = up_i
        self.down_i = down_i
        self.up_var_ratio = up_var_ratio
        self.down_var_ratio = down_var_ratio

    
    def extend_stretches(self, arr, min_seed_len=2, max_gap=1, min_stripe_len=5):

        arr.sort()
        pieces = np.split(arr, np.where(np.diff(arr)!=1)[0]+1)
        filtered = [p for p in pieces if len(p) >= min_seed_len] # can't be empty
        stripes = []
        seed = filtered[0]
        for p in filtered[1:]:
            if p[0] - seed[-1] < (max_gap + 2):
                seed = np.r_[seed, p]
            else:
                if seed[-1] - seed[0] + 1 >= min_stripe_len:
                    stripes.append([seed[0], seed[-1]+1])
                seed = p
        
        if seed[-1] - seed[0] + 1 >= min_stripe_len:
            stripes.append([seed[0], seed[-1]+1])
        
        return stripes
    
    def locate(self, pca1, left_most=True):
        
        loci_1 = self.extend_stretches(np.where(pca1>=0)[0])
        loci_2 = self.extend_stretches(np.where(pca1<=0)[0])
        loci = loci_1 + loci_2
        loci.sort()
        if left_most:
            b = loci[0][1] - 1
        else:
            b = loci[-1][0]
        
        return b
    
    def _extract_xy(self, valid_pixel, maxdis=2000000, mink=3):

        x, y = np.where(valid_pixel)
        M = self.fusion_matrix
        maxdis = maxdis // self.res

        dis = y - x
        local_exp = {}
        for i in range(2, M.shape[0]-1):
            if i > maxdis:
                break
            mask = dis == i
            if mask.sum() >= mink:
                local_x, local_y = x[mask], y[mask]
                tmp_arr = M[local_x, local_y]
                tmp = tmp_arr.mean()
                if tmp > 0:
                    local_exp[i] = tmp
        
        return local_exp
    
    
    def linear_regression(self, exp1, exp2, min_samples=5):

        X = []; Y = []; Xi = []
        for i in sorted(exp1):
            if i in exp2:
                Xi.append(i)
                X.append(exp2[i])
                Y.append(exp1[i])
        X = np.r_[X]
        Y = np.r_[Y]
        Xi = np.r_[Xi]

        if X.size < min_samples:
            rscore = 0
            slope = 0
            warning = False
        else:
            # clean the inputs by isotonic regression
            warning, increasing_bool = check_increasing(Xi, Y)
            X_clean, Y_clean = clean_inputs(Xi, X, Y, increasing_bool)
            rscores = []
            slopes = []
            for X_, Y_ in ([X, Y], [X_clean, Y_clean]):
                X_ = X_[:, np.newaxis]
                huber = HuberRegressor().fit(X_, Y_)
                inlier_mask = np.logical_not(huber.outliers_)
                if inlier_mask.sum() < min_samples:
                    rscore = 0
                    slope = 0
                else:
                    sX = X_[inlier_mask]
                    sY = Y_[inlier_mask]
                    rscore = huber.score(sX, sY)
                    slope = huber.coef_[0]
                rscores.append(rscore)
                slopes.append(slope)

            rscores = np.r_[rscores]
            slopes = np.r_[slopes]
            rscore = rscores.max()
            slope = slopes[np.argmax(rscores)]

        return rscore, slope, warning


    def allele_slope(self):

        # this method assumes detect_bounds has been called
        self.extract_bias(col=self.balance_type)
        shape = self.fusion_matrix.shape
        Bias_M = self.bias[:,np.newaxis] * self.bias
        valid_bias = Bias_M > 0
        p1L = self.fusion_point

        for up_i, down_i in ((self.up_i, self.down_i), (0, self.bias.size-1)):
             ## test bunch of parameters
            u_g_rscores = []
            u_g_slopes = []
            d_g_rscores = []
            d_g_slopes = []
            m_g_rscores = []
            m_g_slopes = []
            for maxdis in [200000, 400000, 500000, 1000000, 1500000, 2000000]:
                for mink in [3]:
                    # middle
                    inter_mask = np.zeros(shape, dtype=bool)
                    inter_mask[up_i:p1L, p1L:down_i] = True
                    valid_pixel = inter_mask & valid_bias
                    m_exp = self._extract_xy(valid_pixel, maxdis=maxdis, mink=mink)

                    # upstream
                    inter_mask = np.zeros(shape, dtype=bool)
                    inter_mask[up_i:p1L, up_i:p1L] = True
                    valid_pixel = inter_mask & valid_bias
                    u_exp = self._extract_xy(valid_pixel, maxdis=maxdis, mink=mink)

                    # downstream
                    inter_mask = np.zeros(shape, dtype=bool)
                    inter_mask[p1L:down_i, p1L:down_i] = True
                    valid_pixel = inter_mask & valid_bias
                    d_exp = self._extract_xy(valid_pixel, maxdis=maxdis, mink=mink)

                    # upstream vs gm expected
                    r, s, w = self.linear_regression(u_exp, self.expected)
                    u_g_rscores.append(r)
                    u_g_slopes.append(s)

                    # downstream vs gm expected
                    r, s, w = self.linear_regression(d_exp, self.expected)
                    d_g_rscores.append(r)
                    d_g_slopes.append(s)

                    # middle vs gm expected
                    r, s, w = self.linear_regression(m_exp, self.expected)
                    m_g_rscores.append(r)
                    m_g_slopes.append(s)
            
            u_g_rscores = np.r_[u_g_rscores]
            d_g_rscores = np.r_[d_g_rscores]
            m_g_rscores = np.r_[m_g_rscores]
            self.u_g_s = u_g_slopes[u_g_rscores.argmax()]
            self.d_g_s = d_g_slopes[d_g_rscores.argmax()]
            self.m_g_s = m_g_slopes[m_g_rscores.argmax()]
            
            self.u_g_r = u_g_rscores.max()
            self.d_g_r = d_g_rscores.max()
            self.m_g_r = m_g_rscores.max()

            if (self.u_g_r > 0.6) and (self.d_g_r > 0.6) and (self.m_g_r > 0.6) and (self.m_g_s > 0):
                break

    
    def correct_heterozygosity(self):

        ## This method assumes allele slope has been called
        junc = self.fusion_point
        hcm = self.fusion_matrix.copy() # heterozygosity corrected matrix
        log.info('{0}: left slope: {1}, R^2: {2}'.format(self.name, self.u_g_s, self.u_g_r))
        log.info('{0}: right slope: {1}, R^2: {2}'.format(self.name, self.d_g_s, self.d_g_r))
        log.info('{0}, middle slope: {1}, R^2: {2}'.format(self.name, self.m_g_s, self.m_g_r))
        self.u_s = self.d_s = 0
        if (self.u_g_r > 0.6) and (self.d_g_r > 0.6):
            slope = max(self.u_g_s, self.d_g_s)
            ms = min(self.u_g_s, self.d_g_s)
            hcm[:junc, :junc] = hcm[:junc, :junc] / slope
            hcm[junc:, junc:] = hcm[junc:, junc:] / slope
            self.u_s = self.m_g_s / self.u_g_s
            self.d_s = self.m_g_s / self.d_g_s
            if (self.m_g_r > 0.6) and (self.m_g_s / ms > 0.1):
                factor = min(1/self.m_g_s, 5/slope)
                hcm[self.up_i:junc, junc:self.down_i] = hcm[self.up_i:junc, junc:self.down_i] * factor

        self.hcm = hcm
    
    def remove_gaps(self):

        # call this method after correct_heterozygosity
        M = self.hcm
        w = 5
        # output a gap-free triu matrix
        mask = M.sum(axis = 0) == 0
        index = list(np.where(mask)[0])
        convert_index = np.where(np.logical_not(mask))[0]
        dim = convert_index.size # pseudo chromosome length
        if w * 2 + 1 > dim: # if too noisy, use original matrix
            index = []
            convert_index = np.arange(len(M))
        temp = np.delete(M, index, 0)
        nM = np.delete(temp, index, 1)
        forward_map = dict(zip(np.arange(len(nM)), convert_index)) # gap-free to original

        self.gap_free = nM
        self.index_map = forward_map
    
    
class Peakachu():

    def __init__(self, matrix, upper=4000000, res=10000, protocol='insitu'):

        self.w = 7
        upper = upper // res
        R, C = matrix.nonzero()
        data = matrix[R, C]
        validmask = np.isfinite(data) & (C-R > (-2*self.w)) & (C-R < (upper+2*self.w))
        R, C, data = R[validmask], C[validmask], data[validmask]
        self.M = sparse.csr_matrix((data, (R, C)), shape=matrix.shape)
        self.get_candidates(4, upper)
        self.r = res
        self.protocol = protocol
    
    def get_candidates(self, lower, upper):

        exp_obs = []
        indices_obs = []
        for i in range(lower, upper+1):
            diag = self.M.diagonal(i)
            if diag.size > 5:
                exp_obs.append(diag.mean())
                indices_obs.append(i)
        
        IR = IsotonicRegression(increasing=False, out_of_bounds='clip')
        IR.fit(np.r_[indices_obs], np.r_[exp_obs])
        d = np.arange(lower, upper+1)
        exp = IR.predict(d)
        expected = dict(zip(d, exp))

        x_arr = np.array([], dtype=int)
        y_arr = np.array([], dtype=int)
        idx = np.arange(self.M.shape[0])
        for i in range(lower, upper+1):
            diag = self.M.diagonal(i)
            e = expected[i]
            if (diag.size > 1) and (e > 0):
                xi = idx[:-i]
                yi = idx[i:]
                diag = diag / e
                mask = diag > 0
                x_arr = np.r_[x_arr, xi[mask]]
                y_arr = np.r_[y_arr, yi[mask]]
        
        self.ridx, self.cidx = x_arr, y_arr
    
    def load_models(self, model_fil_path):
        
        import joblib
        
        self.model = joblib.load(model_fil_path)
        self.w = int((np.sqrt(self.model.feature_importances_.size) - 1) / 2)
    
    def load_expected(self, expected_values):

        self.exp_arr = np.r_[[expected_values[i] for i in sorted(expected_values)]]
    
    def getwindow(self, coords, w):
        '''
        Extract features from center pixels.
        '''
        if len(coords) < 2:
            return [], []
            
        coords = np.r_[coords]
        xi, yi = coords[:,0], coords[:,1]
        mask = (xi - w >= 0) & (yi + w + 1 <= self.M.shape[0]) & (yi - xi - 2 >= w)
        xi, yi = xi[mask], yi[mask]
        if xi.size < 2:
            return [], []
        seed = np.arange(-w, w+1)
        delta = np.tile(seed, (seed.size, 1))
        xxx = xi.reshape((xi.size, 1, 1)) + delta.T
        yyy = yi.reshape((yi.size, 1, 1)) + delta
        v = np.array(self.M[xxx.ravel(), yyy.ravel()]).ravel()
        vvv = v.reshape((xi.size, seed.size, seed.size))
        windows, clist = distance_normalize(vvv, self.exp_arr, xi, yi, w)
        fea = []
        for arr in windows:
            tmp = gaussian_filter(arr, sigma=1, order=0)
            scaled_arr = image_normalize(tmp)
            fea.append(scaled_arr.ravel())

        fea = np.r_[fea]
        clist = np.r_[clist]
        
        return fea, clist
    
    def predict(self, thre, no_pool=False, min_count=2, index_map=None, chains=None):
        
        coords = [(r, c) for r, c in zip(self.ridx, self.cidx)]
        loop_list = set()
        fea, clist = self.getwindow(coords, self.w)
        if len(fea) < 2:
            return list(loop_list)
        
        probas = self.model.predict_proba(fea)[:, 1]
        pM = sparse.csr_matrix((probas, (clist[:,0], clist[:,1])),
                            shape=self.M.shape)
        pfilter = probas >= thre
        ri = clist[:, 0][pfilter]
        ci = clist[:, 1][pfilter]
        Donuts = {(i, j): self.M[i,j] for i, j in zip(ri, ci)}
        if not no_pool:
            tmp = self.local_clustering(Donuts, min_count=min_count, index_map=index_map,
                                        chains=chains)
        else:
            tmp = Donuts
        for i, j in tmp:
            loop_list.add((i, j, pM[i,j]))
        
        loop_list = list(loop_list)

        return loop_list
    
    def _local_search(self, candidates, Donuts, r=10000, index_map=None,
        chains=None):

        final_list = set()
        r = max(1, r//self.r)
        for i, j in candidates:
            if (i - r < 0) or (j + r + 1 > self.M.shape[0]):
                final_list.add((i, j))
            else:
                sub = self.M[i-r:i+r+1, j-r:j+r+1].toarray()
                maxi, maxj, v = i, j, Donuts[(i, j)]
                x, y = np.where(sub > v)
                for xi, yi in zip(x, y):
                    mi, mj = i - (r - xi), j - (r - yi)
                    if index_map is None:
                        if (sub[xi, yi] > v) and (chains[mi]==chains[i]) and (chains[mj]==chains[j]):
                            v = sub[xi, yi]
                            maxi, maxj = mi, mj
                    else:
                        if (sub[xi, yi] > v) and (chains[index_map[mi]]==chains[index_map[i]]) and \
                            (chains[index_map[mj]]==chains[index_map[j]]):
                            v = sub[xi, yi]
                            maxi, maxj = mi, mj
                final_list.add((maxi, maxj))
        
        return final_list
    
    def local_clustering(self, Donuts, min_count=2, r=15000, index_map=None,
        chains=None):

        byregion = {}
        if chains is None:
            byregion[(0, 0)] = Donuts
        else:
            for i, j in Donuts:
                if index_map is None:
                    x, y = i, j
                else:
                    x, y = index_map[i], index_map[j]
                c1, c2 = chains[x], chains[y]
                if not (c1, c2) in byregion:
                    byregion[(c1, c2)] = {}
                byregion[(c1, c2)][(i, j)] = Donuts[(i, j)]
        
        coords = set()
        for k in byregion:
            Donuts = byregion[k]
            tmp = self._local_clustering(Donuts, min_count=min_count, r=r)
            '''
            precise_ = self._local_search(tmp, Donuts, r=r, index_map=index_map,
                                        chains=chains)
            '''
            precise_ = set(tmp)
            coords.update(precise_)
        
        return coords


    def _local_clustering(self, Donuts, min_count=1, r=15000):

        final_list = []
        x = np.r_[[i[0] for i in Donuts]]
        y = np.r_[[i[1] for i in Donuts]]
        if x.size == 0:
            return final_list

        x_anchors = self.find_anchors(x, min_count=min_count, min_dis=r)
        y_anchors = self.find_anchors(y, min_count=min_count, min_dis=r)
        r = max(r//self.r, 2)
        visited = set()
        lookup = set(zip(x, y))
        for x_a in x_anchors:
            for y_a in y_anchors:
                sort_list = []
                for i in range(x_a[1], x_a[2]+1):
                    for j in range(y_a[1], y_a[2]+1):
                        if (i, j) in lookup:
                            sort_list.append((Donuts[(i,j)], (i,j)))
                sort_list.sort(reverse=True)
                self._cluster_core(sort_list, r, visited, final_list)
        
        sort_list = [] # out of anchor
        for i, j in zip(x, y):
            if (i,j) in visited:
                continue
            sort_list.append((Donuts[(i,j)], (i,j)))
        sort_list.sort(reverse=True)
        self._cluster_core(sort_list, r, visited, final_list)
        
        for i, j in zip(x, y):
            if (i,j) in visited:
                continue

            final_list.append((int(i), int(j)))
        
        final_list = list(set(final_list))
        
        return final_list
    
    def find_anchors(self, pos, min_count=3, min_dis=15000, wlen=40000):

        from collections import Counter
        from scipy.signal import find_peaks, peak_widths

        min_dis = max(min_dis//self.r, 2)
        if self.r <= 10000:
            wlen = min(wlen//self.r, 20)
        else:
            wlen = 4

        count = Counter(pos)
        refidx = range(min(count), max(count)+1)
        signal = np.r_[[count[i] for i in refidx]]
        summits = find_peaks(signal, height=min_count, distance=min_dis)[0]
        sorted_summits = [(signal[i],i) for i in summits]
        sorted_summits.sort(reverse=True) # sort by peak count
        
        peaks = set()
        records = {}
        for _, i in sorted_summits:
            tmp = peak_widths(signal, [i], rel_height=1, wlen=wlen)[2:4]
            li, ri = int(np.round(tmp[0][0])), int(np.round(tmp[1][0]))
            lb = refidx[li]
            rb = refidx[ri]
            if not len(peaks):
                peaks.add((refidx[i], lb, rb))
                for b in range(lb, rb+1):
                    records[b] = (refidx[i], lb, rb)
            else:
                for b in range(lb, rb+1):
                    if b in records:
                        # merge anchors
                        m_lb = min(lb, records[b][1])
                        m_rb = max(rb, records[b][2])
                        summit = records[b][0] # always the highest summit
                        peaks.remove(records[b])
                        break
                else: # loop terminates normally
                    m_lb, m_rb, summit = lb, rb, refidx[i]
                peaks.add((summit, m_lb, m_rb))
                for b in range(m_lb, m_rb+1):
                    records[b] = (summit, m_lb, m_rb)
        
        return peaks

    def _cluster_core(self, sort_list, r, visited, final_list):

        from sklearn.cluster import dbscan
        from scipy.spatial.distance import euclidean

        pos = np.r_[[i[1] for i in sort_list]]
        if len(pos) >= 2:
            _, labels = dbscan(pos, eps=r, min_samples=2)
            pool = set()
            for i, p in enumerate(sort_list):
                if p[1] in pool:
                    continue
                c = labels[i]
                if c==-1:
                    continue
                sub = pos[labels==c]
                cen = p[1]
                rad = r
                Local = [p[1]]
                ini = -1
                while len(sub):
                    out = []
                    for q in sub:
                        if tuple(q) in pool:
                            continue
                        tmp = euclidean(q, cen)
                        if tmp<=rad:
                            Local.append(tuple(q))
                        else:
                            out.append(tuple(q))
                    if len(out)==ini:
                        break
                    ini = len(out)
                    tmp = np.r_[Local]
                    # assign centroid to a certain pixel
                    cen = tuple(tmp.mean(axis=0).round().astype(int))
                    rad = np.int64(np.round(max([euclidean(cen,q) for q in Local]))) + r
                    sub = np.r_[out]
                for q in Local:
                    pool.add(q)
                final_list.append((int(p[1][0]), int(p[1][1])))
            
            visited.update(pool)
        
        elif len(pos)==1:
            final_list.append((int(pos[0][0]), int(pos[0][1])))


def combine_annotations(byres):
    """
    Combine loops at different resolutions.

    """
    reslist = sorted(byres)
    pool = byres[reslist[0]]
    visited = [reslist[0]]

    for res in reslist[1:]:
        cur = byres[res]
        for l in cur:
            for pr in visited:
                count = 0
                for i in range(l[1]//pr, l[2]//pr+1):
                    for j in range(l[4]//pr, l[5]//pr+1):
                        tmp = (l[0], i*pr, i*pr+pr, l[3], j*pr, j*pr+pr)
                        if tmp in pool:
                            count += 1
                            pool[tmp].extend(cur[l])
                if count == 0:
                    pool[l] = cur[l]
                     
        visited.append(res)
    
    loop_list = []
    for k in sorted(pool):
        trace = set(pool[k])
        labels = []
        for t in trace:
            labels.append(','.join([t[0], str(t[1]), str(t[2])]))
        label = ','.join(labels)
        loop_list.append(list(map(str, k))+[label])
    
    return loop_list