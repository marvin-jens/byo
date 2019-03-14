from sequence_data.track import Accessor
from logging import debug,warning,error
from sequence_data.io import fasta_chunks
import numpy


def msf_chunks(src):
    """loads a multi-species fasta with one line for each species. By default the first species is the reference and demarcates a new block"""
    ref = ""
    
    species = []
    seqs = []
    extra = ""
    
    for line in src:
        if not ref and line.startswith(">"):
            ref = line.split('.')[0]
            #print "detected ref",ref

        if ref and line.startswith(ref):
            if species and seqs:
                yield zip(species,seqs),extra
                seqs = []
            species.append(ref[1:])
            extra = line.split(ref)[1].rstrip()
        else:
            if line.startswith('>'):
                species.append(line[1:].rstrip())
            else:
                seqs.append(line.rstrip())
    
    if species and seqs:
        yield zip(species,seqs),extra
    
class MSFaAccessor(Accessor):
    """
    multispecies fasta accessor. Computes conservation scores for each column 
    in a whole genome multi-species fasta file (as comes out of GALAXY after 
    stitching MAF blocks)
    """

    def __init__(self,path,chrom,sense,**kwargs):
        debug("# MSFaAccessor: Loading genomic sequence for chromosome %s from '%s'" % (chrom,path))
        fname = os.path.join(path,chrom+".stitched")
        self.chrom = chrom
        try:
            self.data = []
            self.species = []
            
            for fa_id,seq in fasta_chunks(file(fname)):
                self.species.append(fa_id)
                self.data.append(seq.upper())
                debug("# MSFaAccessor: loaded %d bases for '%s'" % (len(seq),fa_id))

            self.n_species = len(self.species) - 1
        except IOError:
            warning("Could not access '%s'. Switching to dummy mode (only zeros)" % fname)
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy

        # TODO: maybe remove this if not needed
        self.get = self.get_oriented
        self.min_cons = numpy.nan
        
    def phylo_score(self,start,end):
        score = numpy.zeros(end-start) + self.min_cons
        for x,col in enumerate(zip(*[s[start:end] for s in self.data])):

            sym = Counter(col[1:])
            
            gaps = sym['-'] 
            match = sym[col[0]]
            mismatch = self.n_species - gaps - match
            
            if gaps < self.n_species:
                score[x] = match - mismatch

        return score
        
    def get_data(self,start,end,sense):
        if start < 0 or end < 0:
            return self.get_dummy(start,end,sense)
        #UCSC convention: start with 1, end is inclusive
        return self.phylo_score(start,end)
        
    def get_oriented(self,start,end,sense):
        if end < 0:
            return self.get_dummy(start,end,sense)
        elif start < 0:
            return numpy.concatenate(self.get_dummy(start,0,sense) + self.get_oriented(0,end,sense))
        if sense == "+":
            return self.phylo_score(start,end)
        else:
            return self.phylo_score(start,end)[::-1]

    def get_dummy(self,start,end,sense):
        return numpy.zeros(end-start) + self.min_cons
