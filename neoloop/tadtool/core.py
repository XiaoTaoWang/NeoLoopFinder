#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import logging
from tadlib.domaincaller.chromLev import Chrom as minChrom
from scipy.sparse import csr_matrix
import numpy as np

log = logging.getLogger(__name__)

    
class TADcaller(minChrom):

    def __init__(self, hicdata, res, model, window_size=2000000):

        hicdata = csr_matrix(hicdata, shape=hicdata.shape)
        minChrom.__init__(self, 'chrN', res, hicdata)
        ws = min(window_size // res, hicdata.shape[0]-5)
        self._dw = ws
        self.windows = np.ones(hicdata.shape[0], np.int32) * ws
        self.hmm = model
    
    def callDomains(self):

        self.calDI(self.windows, 0)
        self.splitChrom(self.DIs)
        minDomains = self.minCore(self.regionDIs)
        self.domains = self.getDomainList(minDomains)
    
    def loop_like(self):

        L = []
        for d in self.domains:
            L.append((d[0]//self.res, d[1]//self.res, 0))
        
        L.sort()

        return L

