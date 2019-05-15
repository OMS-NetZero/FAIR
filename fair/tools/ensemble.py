from __future__ import division

import numpy as np
import scipy.stats as st
import os
from functools import reduce

def mvlognorm(data, n=1000, seed=None, correlated=True):
    """Returns joint lognormal random variables.
    
    Inputs:
        data: (m, p) array of 'observations' where m is the number
            of observations and p is the number of parameters to
            estimate.
        
    Keywords:
        n: number of 
        seed: random seed for variable generation
        correlated: logical. If True, assume random variables are
            correlated, and calculate the correlation coefficient.
        
    Outputs:
        (n, p) array of simulated joint lognormal random variables.
        
    Method: This function takes the input data (untransformed) and 
    calculates the mean, standard deviation and correlation coefficient
    of the data, transforms it to a lognormal, and returns simulated
    data based on the input 'observations'.
    
    It is based on the MethylCapSig R package by Deepak N. Ayyala et al.,
    https://CRAN.R-project.org/package=MethylCapSig
    """
    
    p = data.shape[1]
    mu = np.mean(data, axis=0)
    Sigma = np.var(data, axis=0)
    if correlated:
        corr = np.corrcoef(np.log(data), rowvar=False)
    else:
        corr = np.eye((p))
    
    alpha = np.log(mu) - 0.5*np.log(1.0 + Sigma/mu**2)
    beta = np.diag(np.log(1.0 + Sigma/(mu**2)))
    
    delta = reduce(np.matmul, [np.sqrt(beta), corr, np.sqrt(beta)])
    delta_eigenvalues, delta_eigenvectors = np.linalg.eig(delta)
    root_delta = reduce(np.matmul, 
      [delta_eigenvectors, np.diag(np.sqrt(delta_eigenvalues)), delta_eigenvectors.T])

    out = st.norm.rvs(size=(n,p), random_state=seed)
    for i in range(n):
        out[i,:] = np.exp(alpha[:] + np.matmul(root_delta, out[i,:]))

    return out


def tcrecs_generate(tcrecs_in='cmip5', dist='lognorm', n=1000, correlated=True,
                    strip_ecs_lt_tcr=True,
                    seed=None):
    """Generates a distribution of TCR and ECS.
    
    Inputs:
        tcrecs_in: either 'cmip5' for pre-shipped CMIP5 TCR (transient climate
        response) and ECS (equilibrium climate sensitivity) values, or a 
        2-column array of TCR and ECS values to sample from.
        
    Keywords:
        dist: string. Distribution to use when constructing the
            joint distribution. Accepted values are norm and 
            lognorm (default).
        n: number of samples to generate. Default 1000.
        correlated: logical. If True (default), assume ECS and TCR
            inputs are correlated. The function calculates the
            correlation coefficient automatically.
        strip_ecs_lt_tcr: logical. If True (default), remove values
            where ECS < TCR in place (but still return n samples.)
        seed: random seed for generating variables.

    Output:
        (n, 2) array of sampled ECS, TCR pairs."""
    
    if type(tcrecs_in) is str and tcrecs_in=='cmip5':
        filepath = os.path.join(os.path.dirname(__file__),
          'tcrecs/cmip5tcrecs.csv')
        tcrecs_in = np.loadtxt(filepath, delimiter=',', skiprows=3)
    try:
        assert(type(tcrecs_in) is np.ndarray)
        assert(tcrecs_in.ndim == 2)
        assert(tcrecs_in.shape[1] == 2)
    except AssertionError:
        raise ValueError('tcrecs_in should "cmip5" or an array of shape (n, 2)')
    
    dist = dist.lower()
    
    def _genvar(tcrecs, dist, n, seed, correlated):
        if dist=='lognorm':
            out = mvlognorm(tcrecs, n=n, seed=seed, correlated=correlated)
        elif dist=='norm':
            mu = np.mean(tcrecs, axis=0)
            if correlated:
                cov = np.cov(tcrecs, rowvar=False)
            else:
                cov = np.diag(np.var(tcrecs, axis=0))
            out = st.multivariate_normal.rvs(mu, cov, size=n, random_state=seed)
        else:
            raise ValueError('dist should be "norm" or "lognorm"')
        return out

    tcrecs_out = _genvar(tcrecs_in, dist, n, seed, correlated)
    
    if strip_ecs_lt_tcr:
        tcrecs_out = np.delete(tcrecs_out, 
          np.where(tcrecs_out[:,0] > tcrecs_out[:,1]), axis=0)
        nreq = n - len(tcrecs_out[:,0])
        while nreq>0:
            # required otherwise we are repeating values
            # this is still deterministic for constant ensemble size
            if seed is not None:
                seed = seed + nreq
            new = _genvar(tcrecs_out, dist, n-len(tcrecs_out[:,0]),
              seed, correlated)
            tcrecs_out = np.append(tcrecs_out, new, axis=0)
            tcrecs_out = np.delete(tcrecs_out, np.where(tcrecs_out[:,0] > 
              tcrecs_out[:,1]), axis=0)
            nreq = n - len(tcrecs_out[:,0])
    return tcrecs_out
