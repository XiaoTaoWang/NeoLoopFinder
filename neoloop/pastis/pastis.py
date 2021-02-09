#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import string
import numpy as np
from neoloop.assembly import complexSV
from scipy import optimize, sparse
from neoloop.pastis import possion_model

def MDS_obj(X, distances):

    X = X.reshape(-1, 3)
    dis = np.sqrt(((X[distances.row] - X[distances.col])**2).sum(axis=1))

    return ((dis - distances.data)**2 / distances.data**2).sum()

def MDS_gradient(X, distances):

    X = X.reshape(-1, 3)
    dis = np.sqrt(((X[distances.row] - X[distances.col])**2).sum(axis=1))

    grad = 2 * ((dis - distances.data) / dis /
                distances.data**2)[:, np.newaxis] * (
        X[distances.row] - X[distances.col])
    grad_ = np.zeros(X.shape)

    for i in range(X.shape[0]):
        grad_[i] += grad[distances.row == i].sum(axis=0)
        grad_[i] -= grad[distances.col == i].sum(axis=0)

    X = X.flatten()

    return grad_.flatten()

def poisson_obj(X, counts, alpha=-3., beta=1., bias=None):

    if bias is None:
        bias = np.ones(counts.shape[0])
    bias = bias.flatten()
    
    dis = np.sqrt(((X[counts.row] - X[counts.col])**2).sum(axis=1))
    fdis = bias[counts.row] * bias[counts.col] * beta * dis ** alpha

    obj = fdis.sum() - (counts.data * np.log(fdis)).sum()
    if np.isnan(obj):
        raise ValueError("Objective function is nan")

    return obj

def poisson_gradient(X, counts, alpha=-3., beta=1., bias=None):

    if bias is None:
        bias = np.ones(counts.shape[0])
    bias = bias.flatten()

    dis = np.sqrt(((X[counts.row] - X[counts.col])**2).sum(axis=1))
    fdis = bias[counts.row] * bias[counts.col] * beta * dis ** alpha

    diff = X[counts.row] - X[counts.col]

    grad = - ((counts.data / fdis - 1) * fdis * alpha /
              (dis ** 2))[:, np.newaxis] * diff

    grad_ = np.zeros(X.shape)

    for i in range(X.shape[0]):
        grad_[i] += grad[counts.row == i].sum(axis=0)
        grad_[i] -= grad[counts.col == i].sum(axis=0)

    return grad_


def eval_f(x, user_data=None):

    n, counts, alpha, beta, bias = user_data
    x = x.reshape((n, 3))
    obj = poisson_obj(x, counts, alpha=alpha, beta=beta, bias=bias)
    x = x.flatten()

    return obj


def eval_grad_f(x, user_data=None):

    n, counts, alpha, beta, bias = user_data
    x = x.reshape((n, 3))
    grad = poisson_gradient(x, counts, alpha=alpha, beta=beta, bias=bias)
    x = x.flatten()

    return grad.flatten()


class PMCore(complexSV):

    chain_identifiers = dict(zip(range(26), string.ascii_uppercase))

    def __init__(self, clr, candidate, protocol='insitu', span=None, slopes={}, col='sweight',
        alpha=-3., beta=1., max_iter=10000, max_iter_outer_loop=5,
        seed=0):

        # reconstruct matrix
        if not col in ['weight', 'sweight']:
            if not span is None:
                wk = complexSV(clr, candidate, span=span, flexible=False, col=False,
                               protocol=protocol, slopes=slopes)
            else:
                wk = complexSV(clr, candidate, col=False, protocol=protocol, slopes=slopes)
            wk.reorganize_matrix()
            counts = wk.fusion_matrix
            x, y = counts.nonzero()
            counts[y, x] = counts[x, y]
        else:
            if not span is None:
                wk = complexSV(clr, candidate, span=span, flexible=False, col=col, protocol=protocol,
                               slopes=slopes)
                _wk = complexSV(clr, candidate, span=span, flexible=False, col=False, protocol=protocol,
                                slopes=slopes)
            else:
                wk = complexSV(clr, candidate, col=col, protocol=protocol, slopes=slopes)
                _wk = complexSV(clr, candidate, col=False, protocol=protocol, slopes=slopes)
            wk.reorganize_matrix()
            wk.correct_heterozygosity()
            counts = wk.hcm

            _wk.reorganize_matrix()
            raw = _wk.fusion_matrix # this is an upper triangular matrix
            marg = raw.sum(0) + raw.sum(1) - raw.diagonal(0)
            scale = marg[marg > 0].mean()

            counts = counts * scale
        
        counts = sparse.coo_matrix(counts)
        counts.setdiag(0)
        counts.eliminate_zeros()
        self.fusion_matrix = counts

        # PM2 parameters
        self.alpha = alpha
        self.beta = beta
        self.max_iter = max_iter
        self.max_iter_outer_loop = max_iter_outer_loop
        self.random_state = np.random.RandomState(seed=seed)

        # mark different chains (when write a PDB file)
        self.chains = {}
        for i, r in enumerate(wk.index):
            for ri in range(r[0], r[1]):
                self.chains[ri] = self.chain_identifiers[i]

    
    def _compute_wish_distances(self):

        if self.beta == 0:
            raise ValueError("beta cannot be equal to 0.")

        wish_distances = self.fusion_matrix / self.beta
        wish_distances.data[wish_distances.data != 0] **= 1. / self.alpha

        return wish_distances

    
    def _estimate_X(self):

        n = self.fusion_matrix.shape[0]
        ini = 1 - 2 * self.random_state.rand(n * 3)
        distances = self._compute_wish_distances()
        results = optimize.fmin_l_bfgs_b(
            MDS_obj, ini.flatten(),
            MDS_gradient,
            (distances,),
            iprint=0,
            factr=1e12,
            maxiter=self.max_iter
        )

        return results[0].reshape(-1, 3)

    
    def estimate_X(self, ini):

        n = self.fusion_matrix.shape[0]

        ini = np.array(ini)

        data = (n, self.fusion_matrix, self.alpha_, self.beta_, None)

        results = optimize.fmin_l_bfgs_b(
            eval_f,
            ini.flatten(),
            eval_grad_f,
            (data, ),
            iprint=0,
            maxiter=self.max_iter,
            )

        results = results[0].reshape(-1, 3)

        return results
    
    def fit(self):

        X = self._estimate_X() # initialization strategy: MDS2
        self.alpha_ = self.alpha
        self.beta_ = self.beta
        for _ in range(self.max_iter_outer_loop):
            self.alpha_, self.beta_ = possion_model.estimate_alpha_beta(
                self.fusion_matrix,
                X, bias=None, ini=[self.alpha_, self.beta_],
                random_state=self.random_state
            )
            X_ = self.estimate_X(ini=X)
        
        self.coords = X_
    
    def output_pdb(self, filename, ref_dis=1.1):

        from scipy.spatial.distance import euclidean

        X = self.coords
        q_dis = np.mean([euclidean(X[i+1], X[i]) for i in range(len(X)-1)])
        X = X * ref_dis / q_dis
        with open(filename, 'w') as out:
            for idx in range(X.shape[0]):
                # 1 - 6, Record name
                line = '{0:6s}'.format('ATOM')
                # 7 - 11, Atom serial number
                line += '{0:>5d}'.format(idx+1)
                line += ' '
                # 13 - 16, Atom name
                line += '{0:<4s}'.format('CA')
                # 17, alternate location indicator
                line += ' '
                # 18 - 20, Residue name
                line += 'GLY'
                line += ' '
                # 22, chain identifier
                line += '{0}'.format(self.chains[idx])
                # 23 - 26, residue sequence number
                line += '{0:>4d}'.format(idx+1)
                # 27, insertion code
                line += ' '
                # 28 - 30, blank
                line += '   '
                # 31 - 38
                line += '{0:>8.3f}'.format(X[idx][0])
                # 39 - 46
                line += '{0:>8.3f}'.format(X[idx][1])
                # 47 - 54
                line += '{0:>8.3f}'.format(X[idx][2])
                # 55 - 60, Occupancy
                line += '{0:>6.2f}'.format(1.0)
                # 61 - 66, Temperature factor
                line += '{0:>6.2f}'.format(75.0)
                # 67 - 78, blank
                line += '            '
                # 79 - 80, charge on the atom
                line += ' C\n'
                out.write(line)
            # CONECT section
            for i in range(1, X.shape[0]):
                # 1 - 6, Record name
                line = '{0:6s}'.format('CONECT')
                # 7 - 11, atom serial number
                line += '{0:>5d}'.format(i)
                # 12 - 16, serial number of bonded atom
                line += '{0:>5d}\n'.format(i+1)
                out.write(line)