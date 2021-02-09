#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import numpy as np
from scipy import optimize, sparse

def poisson_exp(X, counts, alpha, bias=None,
                beta=None, use_empty_entries=True):
    
    if bias is None:
        bias = np.ones((counts.shape[0], 1))

    ll = _poisson_exp_sparse(X, counts, alpha, bias=bias, beta=beta,
                            use_empty_entries=use_empty_entries)

    return ll

def _poisson_exp_sparse(X, counts, alpha, bias,
                        beta=None, use_empty_entries=False):
    m, n = X.shape

    d = np.sqrt(((X[counts.row] - X[counts.col])**2).sum(axis=1))
    if use_empty_entries:
        raise NotImplementedError

    bias = bias.flatten()
    if beta is None:
        beta = counts.sum() / (
            (d ** alpha) * bias[counts.row] * bias[counts.col]).sum()

    g = beta * d ** alpha * \
        bias[counts.row] * bias[counts.col]
    ll = (counts.data * np.log(beta) + alpha * counts.data * np.log(d) +
          counts.data * np.log(bias[counts.row] * bias[counts.col]))
    ll -= g

    if np.isnan(ll.sum()):
        raise ValueError("Objective function is Not a Number")

    return -ll.sum()

def gradient_poisson_exp(X, counts, alpha, bias=None,
                         beta=None, use_empty_entries=True):

    if bias is None:
        bias = np.ones((counts.shape[0], 1))

    return _gradient_poisson_exp_sparse(
            X, counts, alpha, bias, beta,
            use_empty_entries=use_empty_entries)

def _gradient_poisson_exp_sparse(X, counts, alpha, bias, beta,
                                 use_empty_entries=True):
    m, n = X.shape
    bias = bias.flatten()

    if use_empty_entries:
        raise NotImplementedError

    d = np.sqrt(((X[counts.row] - X[counts.col])**2).sum(axis=1))

    beta = counts.sum() / (
        (d ** alpha) * bias[counts.row] * bias[counts.col]).sum()

    grad_alpha = - beta * (bias[counts.row] * bias[counts.col] * d ** alpha *
                           np.log(d)).sum() \
        + (counts.data *
           np.log(d)).sum()

    return - np.array([grad_alpha])


def eval_f(x, user_data=None):

    m, n, counts, X, bias, use_empty_entries = user_data
    X = X.reshape((m, n))
    tmp = poisson_exp(X, counts, x[0], bias=bias,
                      use_empty_entries=use_empty_entries)
    return tmp


def eval_grad_f(x, user_data=None):
    """
    Evaluate the gradient of the function in alpha
    """

    m, n, counts, X, bias, use_empty_entries = user_data
    X = X.reshape((m, n))
    tmp = gradient_poisson_exp(X, counts, x[0], bias=bias, beta=None,
                               use_empty_entries=use_empty_entries)

    return tmp

def _estimate_beta(counts, X, alpha=-3, bias=None):

    if bias is None:
        bias = np.ones((counts.shape[0], 1))
    
    dis = np.sqrt(((X[counts.row] - X[counts.col])**2).sum(axis=1))
    bias = bias.flatten()
    beta = counts.sum() / (
            (dis ** alpha) * bias[counts.row] * bias[counts.col]).sum()
    
    return beta


def estimate_alpha_beta(counts, X, bias=None, ini=None,
                        use_empty_entries=False, random_state=None):

    m, n = X.shape
    bounds = np.array(
        [[-100, 1e-2]])

    if ini is None:
        ini = - random_state.randint(1, 100, size=(2, )) + \
            random_state.rand(1)

    data = (m, n, counts, X, bias,
            use_empty_entries)

    results = optimize.fmin_l_bfgs_b(
        eval_f,
        ini[0],
        eval_grad_f,
        (data, ),
        bounds=bounds,
        iprint=1,
        maxiter=1000,
        )

    beta = _estimate_beta(counts, X, alpha=results[0], bias=bias)

    return results[0], beta