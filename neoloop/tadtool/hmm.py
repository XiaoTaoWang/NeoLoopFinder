#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import cooler, logging
from scipy.sparse import triu
from tadlib.hitad.chromLev import Chrom
from pomegranate import NormalDistribution, HiddenMarkovModel, GeneralMixtureModel, State

log = logging.getLogger(__name__)

class Genome(object):

    def __init__(self, uri, balance_type='sweight', cpu_core=1):

        if not balance_type in ['weight', 'sweight']:
            correct = False
        else:
            correct = balance_type
        
        lib = cooler.Cooler(uri)
        res = lib.binsize
        seqs = []
        log.debug('Calculating DIs for each chromosome ...')
        for c in lib.chromnames:
            if not c in ['chr1', 'chr2', 'chr3']:
                continue
            log.debug('Chrom {0} ...'.format(c))
            tdata = triu(lib.matrix(balance=correct, sparse=True).fetch(c)).tocsr()
            work = Chrom(c, res, tdata, 'pseudo', maxapart=3000000)
            work.minWindows(0, work.chromLen, work._dw)
            work.calDI(work.windows, 0)
            work.splitChrom(work.DIs)
            for r in work.regionDIs:
                withzeros = work.regionDIs[r]
                nozeros = withzeros[withzeros!=0]
                if nozeros.size > 20:
                    seqs.append(nozeros)
        
        self.training_data = seqs
        log.debug('Initialize a Hidden Markov Model ...')
        self.model = self.oriHMMParams()
        log.debug('Learning parameters ...')
        self.fit(cpu_core=cpu_core)
    
    def oriHMMParams(self):
        """
        Set initial parameters for the Hidden Markov Model (HMM).
        
        Attributes
        ----------
        HMMParams : dict
            Has 3 keys: "A", state transition matrix, "B" (emission probabilities),
            specifying parameters (Means, Variances, Weights) of the mixture
            Gaussian distributions for each hidden state, and "pi", indicating
            the hidden state weights. This dict will be updated after learning
            procedure.
        """
        hmm = HiddenMarkovModel()
        # GMM emissions
        # 4 Hidden States:
        # 0--start, 1--downstream, 2--upstream, 3--end
        numdists = 3 # Three-distribution Gaussian Mixtures
        var = 7.5 / (numdists - 1)
        means = [[], [], [], []]
        for i in range(numdists):
            means[3].append(i * 7.5 / ( numdists - 1 ) + 2.5)
            means[2].append(i * 7.5 / ( numdists - 1 ))
            means[1].append(-i * 7.5 / ( numdists - 1 ))
            means[0].append(-i * 7.5 / ( numdists - 1 ) - 2.5)
        states = []
        for i, m in enumerate(means):
            tmp = []
            for j in m:
                tmp.append(NormalDistribution(j, var))
            mixture = GeneralMixtureModel(tmp)
            states.append(State(mixture, name=str(i)))
        hmm.add_states(*tuple(states))

        # Transmission matrix
        #A = [[0., 1., 0., 0.],
        #    [0., 0.5, 0.5, 0.],
        #    [0., 0., 0.5, 0.5],
        #    [1., 0., 0., 0.]]
        hmm.add_transition(states[0], states[1], 1)
        hmm.add_transition(states[1], states[1], 0.5)
        hmm.add_transition(states[1], states[2], 0.5)
        hmm.add_transition(states[2], states[2], 0.5)
        hmm.add_transition(states[2], states[3], 0.5)
        hmm.add_transition(states[3], states[0], 1)

        #pi = [0.2, 0.3, 0.3, 0.2]
        hmm.add_transition(hmm.start, states[0], 1)
        hmm.add_transition(states[3], hmm.end, 1)

        hmm.bake()

        return hmm

    
    def fit(self, cpu_core):

        self.model.fit(self.training_data, algorithm='baum-welch', max_iterations=10000,
                  stop_threshold=1e-5, n_jobs=cpu_core, verbose=False)
        

