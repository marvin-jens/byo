import logging
from blist import blist
import os

class SparseMap(object):
    """
    This class implements a memory and time-efficient data structure 
    to store and retrieve features on a linear interval (i.e. a chromosome).
    All the heavy lifting is done by the well-tested `blist` module.
    This class adds load/store for persistence (requires cPickle) and a nice interface.
    """
    
    def __init__(self,L, el_type = set, overwrite=False):
        self.logger = logging.getLogger("SparseMap")
        self.L = L
        self.overwrite = overwrite
        self.logger.info("initializing empty blist with {0} entries".format(L))
        if L:
            self.data = blist([el_type()]) * L
        self.el_type = el_type
        
    def record(self, start, end, item):
        for x in xrange(start, end):
            stored = self.data[x]
            if not stored or self.overwrite:
                # make a new set, because they are 
                # initially all shallow copies of the same set() !
                self.data[x] = self.el_type(item)
            else:
                self.data[x].add(item)

    def __getitem__(self, pos):
        return self.data[pos]
    
    def __getslice__(self, i, j):
        return self.data[i:j]
    
    def __setslice__(self, i, j, value):
        #self.data[i:j] = blist([value]) * (j-i)
        self.record(i, j, value)
        
    def store(self, fname):
        self.logger.info("storing compressed pickle stream to '{0}'".format(fname))
        import cPickle as pickle
        import lz4framed as LZ4
        pick = pickle.dumps(self.data, 2)
        self.logger.debug("pickle done. compressing {0} bytes...".format(len(pick)))
        comp = LZ4.compress(pick)
        self.logger.debug("compression done. writing {0} bytes...".format(len(comp)))
        file(fname,'wb').write(comp)
        self.logger.info("done.")

    def load(self, fname):
        self.logger.info("loading compressed pickle stream from '{0}'".format(fname))
        import cPickle as pickle
        import lz4framed as LZ4
        self.data = pickle.loads(LZ4.decompress(file(fname,'rb').read()))
        self.L = len(self.data)
        self.logger.info("done.")

from byo.track import Accessor
class SparseMapAccessor(Accessor):
    def __init__(self,path,chrom,sense,sense_specific=False,system=None,ext=".lz",**kwargs):
        super(SparseMapAccessor,self).__init__(path,chrom,sense,sense_specific=sense_specific,**kwargs)
        self.logger = logging.getLogger('byo.io.SparseMapAccessor')
        self.logger.debug("loading compressed sparse data for chromosome {0}{1} from '{2}'".format(chrom,sense,path))

        self.sense_specific = sense_specific
        if sense_specific:
            fname = os.path.join(path,chrom+sense+ext)
        else:
            fname = os.path.join(path,chrom+ext)

        self.smap = SparseMap(0)
        self.smap.load(fname)
        
        # link this up directly, so other layers on top of this can treat it like an ArrayAccessor
        self.data = self.smap.data

    def get_data(self, chrom, start, end, sense):
        return self.smap.data[start, end]
        
    
if __name__ == "__main__":
    FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)-25s\t%(message)s'

    logging.basicConfig(level=logging.DEBUG, format=FORMAT)
    import byo.gene_model
    S = SparseMap(249250621)
    S.record(10,100,"a")
    S.record(80,170,"b")
    S.record(110,200,"c")
    S.record(110,140,"d")
    S.record(3000,5000,"e")

    for p in [9,10,99,100,105,300]:
        print p,"->",S[p]
    
    for s,e in [(1,9),(5,15),(100,130)]:
        print s,e, "->",S[s:e]
    
    S.store("s.lz")
    S.load("s.lz")
    for s,e in [(1,9),(5,15),(100,130)]:
        print s,e, "->",S[s:e]
    