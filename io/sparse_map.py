import logging

class SparseMap(object):
    """
    This class implements a memory and time-efficient data structure 
    to store and retrieve features on a linear interval (i.e. a chromosome).
    All the heavy lifting is done by the well-tested `blist` module.
    This class adds load/store for persistence (requires cPickle) and a nice interface.
    """
    
    def __init__(self,L):
        from blist import blist
        self.logger = logging.getLogger("SparseMap")
        self.L = L
        self.logger.info("initializing empty blist with {0} entries".format(L))
        self.data = blist([set()]) * L
        
    def add(self, start, end, item):
        for x in xrange(start, end):
            stored = self.data[x]
            if not stored:
                # make a new set, because they are 
                # initially all shallow copies of the same set() !
                self.data[x] = set(item)
            else:
                self.data[x].add(item)

    def __getitem__(self, pos):
        return self.data[pos]
    
    def __getslice__(self, i, j):
        return self.data[i:j]
    
    def store(self, fname):
        self.logger.info("storing compressed pickle stream to '{0}'".format(fname))
        import cPickle as pickle
        import lz4framed as LZ4
        file(fname,'wb').write(LZ4.compress(pickle.dumps(self.data)))
        self.logger.info("done.")

    def load(self, fname):
        self.logger.info("loading compressed pickle stream from '{0}'".format(fname))
        import cPickle as pickle
        import lz4framed as LZ4
        self.data = pickle.loads(LZ4.decompress(file(fname,'rb').read()))
        self.L = len(self.data)
        self.logger.info("done.")

if __name__ == "__main__":
    FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)-25s\t%(message)s'

    logging.basicConfig(level=logging.DEBUG, format=FORMAT)
    import byo.gene_model
    S = SparseMap(249250621)
    S.add(10,100,"a")
    S.add(80,170,"b")
    S.add(110,200,"c")
    S.add(110,140,"d")
    S.add(3000,5000,"e")

    for p in [9,10,99,100,105,300]:
        print p,"->",S[p]
    
    for s,e in [(1,9),(5,15),(100,130)]:
        print s,e, "->",S[s:e]
    
    S.store("s.lz")
    S.load("s.lz")
    for s,e in [(1,9),(5,15),(100,130)]:
        print s,e, "->",S[s:e]
    