#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

class Loops(object):

    def __init__(self, neoloop_fil):
        
        self.loops = self.parse_loops(neoloop_fil)
    
    def parse_loops(self, filename):

        loops = {}
        with open(filename, 'r') as source:
            for line in source:
                p = line.rstrip().split()
                loci1 = (p[0], int(p[1]), int(p[2]))
                loci2 = (p[3], int(p[4]), int(p[5]))
                label = p[-1].split(',')
                IDs = [label[i] for i in range(0, len(label), 3)]
                neoloop = [int(label[i+2]) for i in range(0, len(label), 3)]
                for i, n in zip(IDs, neoloop):
                    if not i in loops:
                        loops[i] = []
                    loops[i].append((loci1, loci2, n))
        
        return loops
