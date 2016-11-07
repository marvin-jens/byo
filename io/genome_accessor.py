import os
import mmap
import re
import logging
from byo.track import Accessor, Track
from byo import complement, rev_comp

import urllib2

class RemoteCache(object):
    def __init__(self, urlbase):
        self.urlbase = urlbase
        self.logger = logging.getLogger("byo.io.genome_accessor.RemoteCache('{0}')".format(urlbase))

    def __getitem__(self, genome):
        self.logger.debug('getitem {0}'.format(genome))
            
        class Proxy(object):
            no_data = False
            def get_oriented(this, chrom, start, end, strand):
                url = "{base}/{genome}/{chrom}:{start}-{end}{strand}".format(
                    base = self.urlbase, 
                    genome=genome, 
                    chrom=chrom, 
                    start=start, 
                    end=end, 
                    strand=strand
                )
                self.logger.debug('requesting "{0}"'.format(url))
                res = urllib2.urlopen(url)
                return res.read()

        return Proxy()


class Singleton(type):
    """
    To be used as a metaclass. Ensures a class is only instantiated once 
    **per combination of args and kwargs**
    """
    _instances = {}

    def __call__(cls, *args, **kwargs):
        key = (cls, args, tuple(kwargs.items()) )
        if key not in cls._instances:
            cls._instances[key] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[key]




class GenomeCache(object):

    # ensure that only one GenomeCache object is created per path during run-time
    __metaclass__ = Singleton 

    def __init__(self,path):
        self.cached = {}
        self.path = path
        self.logger = logging.getLogger("byo.io.genome_accessor.GenomeCache")

    def __getitem__(self,name):
        if not name in self.cached:
            self.logger.debug("{0} not in genome cache".format(name))
            
            # first look for 2bit file
            self.cached[name] = Track(self.path,TwoBitAccessor,system=name,split_chrom='.')
            if self.cached[name].no_data:
                # not found, look for fasta or LZ compressed 
                self.cached[name] = Track(self.path,GenomeAccessor,system=name,split_chrom='.')
                
            self.logger.info("{0} genomes loaded".format(len(self.cached)))
        
        return self.cached[name]

    def available_genomes(self):
        import glob
        def listit(fmt):
            return [os.path.splitext(
                os.path.basename(f)
            )[0] for f in sorted(glob.glob(self.path+'/*.{0}'.format(fmt)))]
            
        twobit = listit('2bit')
        fasta = listit('fa')
        lz = listit('lz')
        
        return twobit + fasta + lz

class TwoBitAccessor(Accessor):
    def __init__(self, path,chrom,sense,system=None, split_chrom="", **kwargs):
        kwargs['sense_specific'] = False
        super(TwoBitAccessor,self).__init__(path,chrom,sense,system=system,**kwargs)
        import twobitreader as TB
        from time import time

        self.system = system
        self.data = None
        
        # try to access the whole genome, using indexing for fast lookup
        fname = os.path.join(path,"{0}.2bit".format(system))

        self.logger = logging.getLogger('byo.io.TwoBitAccessor({0})'.format(fname))
        t0 = time()
        try:
            self.data = TB.TwoBitFile(fname)
            t1 = time()
            self.data.chrom_stats = self.data.sequence_sizes()
            self.chrom_lookup = {}
            if split_chrom:
                self.chrom_lookup = {}
                for chr in self.data.chrom_stats.keys():
                    short = chr.split(split_chrom)[0]
                    self.chrom_lookup[short] = chr

                self.covered_strands = [chrom+'+' for chrom in self.chrom_lookup.keys()] + [chrom+'-' for chrom in self.chrom_lookup.keys()]
            else:
                self.covered_strands = [chrom+'+' for chrom in self.data.keys()] + [chrom+'-' for chrom in self.data.keys()]
            self.logger.info("file provides {0} sequences.".format(len(self.data.chrom_stats)))

        except IOError:
            t1 = time()
            # all fails: return Ns only
            self.logger.warning("Could not access '{0}'. Switching to dummy mode (only Ns)".format(fname))
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy
            self.covered_strands = [chrom+'+',chrom+'-']
            self.no_data = True

        # TODO: maybe remove this if not needed
        self.get = self.get_oriented
        #self.logger.debug("covered strands: '{0}'".format(",".join(self.covered_strands[:10])) )
        t2 = time()
        self.logger.debug("opening 2bit file took {0:.1f}ms, entire constructor {1:.1f}ms".format( (t1-t0)*1000., (t2-t0)*1000. ) )

    def get_dummy(self, chrom, start, end, sense):
        return "N" *(end-start)

    def get_data(self,chrom,start,end,sense):
        if self.chrom_lookup:
            chrom = self.chrom_lookup[chrom]

        seq = self.data[chrom][start:end]
            
        if sense == "-":
            seq = complement(seq)

        return seq



class IndexedFasta(object):
    def __init__(self,fname,split_chrom="",lz = False, **kwargs):
        self.logger = logging.getLogger('byo.io.IndexedFasta')
        self.fname = fname
        self.chrom_stats = {}
        self.split_chrom = split_chrom

        if lz:
            from byo.io.lz import LZFile
            self.logger.debug("{0} is LZ file".format(fname))
            self._f = LZFile(fname)
            idx_file = self._f
        else:
            f = file(fname)
            idx_file = f
            self._f = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)

        # try to load index
        ipath = fname + '.byo_index'
        if os.access(ipath,os.R_OK):
            self.load_index(ipath)
        else:
            self.index(idx_file)
            self.store_index(ipath)



    def index(self, idx_file):
        self.logger.debug("index('{self.fname}') split_chrom={self.split_chrom}".format(**locals()))

        ofs = 0
        chrom = "undef"
        chrom_ofs = 0
        size = 0
        nl_char = 0

        for line in idx_file:
            ofs += len(line)
            if line.startswith('>'):
                # store size of previous sequence
                if size: self.chrom_stats[chrom].append(size)

                chrom = line[1:].split()[0].strip()
                if self.split_chrom:
                    # use this to strip garbage from chrom/contig name like chr1:new 
                    # ->split_chrom=':' -> chrom=chr1
                    chrom = chrom.split(self.split_chrom)[0]
                chrom_ofs = ofs
            else:
                if not chrom in self.chrom_stats:
                    # this is the first line of the new chrom
                    size = 0
                    lline = len(line)
                    ldata = len(line.strip())
                    nl_char = lline - ldata
                    self.chrom_stats[chrom] = [chrom_ofs,ldata,nl_char,line[ldata:]]
                size += len(line.strip())

        # store size of previous sequence
        if size: self.chrom_stats[chrom].append(size)
        idx_file.flush()
        idx_file.seek(0)

    def store_index(self,ipath):
        self.logger.info("store_index('%s')" % ipath)
        
        # write to tmp-file first and in the end rename in order to have this atomic 
        # otherwise parallel building of the same index may screw it up.
        
        import tempfile
        tmp = tempfile.NamedTemporaryFile(mode="w",dir = os.path.dirname(ipath),delete=False)
        for chrom in sorted(self.chrom_stats.keys()):
            ofs,ldata,skip,skipchar,size = self.chrom_stats[chrom]
            tmp.write("%s\t%d\t%d\t%d\t%r\t%d\n" % (chrom,ofs,ldata,skip,skipchar,size))
        
        # make sure everything is on disk
        os.fsync(tmp)
        tmp.close()
        
        # make it accessible to everyone
        import stat
        os.chmod(tmp.name, stat.S_IROTH | stat.S_IRGRP | stat.S_IRUSR)
        
        # this is atomic on POSIX as we have created tmp in the same directory, 
        # therefore same filesystem
        os.rename(tmp.name,ipath)
        

    def load_index(self,ipath):
        self.logger.info("load_index('%s')" % ipath)
        self.chrom_stats = {}
        for line in file(ipath):
            chrom,ofs,ldata,skip,skipchar,size = line.rstrip().split('\t')
            self.chrom_stats[chrom] = (int(ofs),int(ldata),int(skip),skipchar[1:-1].decode('string_escape'),int(size))
        
    
    def get_data(self,chrom,start,end,sense):
        if not self.chrom_stats:
            self.index()

        ofs,ldata,skip,skip_char,size = self.chrom_stats[chrom]
        #print "ldata",ldata
        #print "chromstats",self.chrom_stats[chrom]
        pad_start = 0
        pad_end = 0
        if start < 0:
            pad_start = -start
            start = 0

        if end > size:
            pad_end = end - size
            end = size

        l_start = start / ldata
        l_end = end / ldata
        #print "lines",l_start,l_end
        ofs_start = l_start * skip + start + ofs
        ofs_end = l_end * skip + end + ofs
        #print "ofs",ofs_start,ofs_end,ofs_end - ofs_start
        
        s = self._f[ofs_start:ofs_end].replace(skip_char,"")
        if pad_start or pad_end:
            s = "N"*pad_start + s + "N"*pad_end

        if sense == "-":
            s = rev_comp(s)
        return s

        
class GenomeAccessor(Accessor):
    def __init__(self,path,chrom,sense,system=None,**kwargs):
        super(GenomeAccessor,self).__init__(path,chrom,sense,system=system,**kwargs)
        self.logger = logging.getLogger('byo.io.GenomeAccessor')
        self.logger.debug("opening genomic sequence for chromosome %s from '%s'" % (chrom,path))

        self.system = system
        self.data = None
        
        # try to access the whole genome, using indexing for fast lookup
        trials = [
            (os.path.join(path,system+'.fa'),'*'),
            (os.path.join(path,system+'.fna'),'*'),
            (os.path.join(path,system+'.fasta'),'*'),
            (os.path.join(path,chrom+".fa"),[chrom+'+',chrom+'-'])
        ]
        for fname,strands in trials:
            self.logger.debug("trying to load '%s'" % fname)
            if os.access(fname + '.lzot', os.R_OK):
                self.logger.debug("{0} is LZ file".format(fname))
                self.data = IndexedFasta(fname,lz = True, **kwargs)
                self.covered_strands = strands
                break
                
            elif os.access(fname,os.R_OK):
                self.data = IndexedFasta(fname,**kwargs)
                self.covered_strands = strands
                break
                
        if not self.data:
            # all fails: return Ns only
            self.logger.warning("Could not access any of '%s'. Switching to dummy mode (only Ns)" % str(trials))
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy
            self.covered_strands = [chrom+'+',chrom+'-']
            self.no_data = True
        #else:
            ## register for all chroms/strands
            #self.covered_strands = [chrom+'+' for chrom in self.data.chrom_stats.keys()] + [chrom+'-' for chrom in self.data.chrom_stats.keys()]

        # TODO: maybe remove this if not needed
        self.get = self.get_oriented


    def load_indexed(self,path):
        ipath = path+'.index'
        if not os.access(ipath,os.R_OK):
            index = self.build_index(path,ipath)
        else:
            index = self.load_index(ipath)

        self.chrom_ofs = index.chrom_ofs
        

    def get_data(self,chrom,start,end,sense):
        seq = self.data.get_data(chrom,start,end,"+")
            
        if sense == "-":
            seq = complement(seq)

        return seq


    def get_dummy(self,chrom,start,end,sense):
        return "N"*int(end-start)



class MSFAccessor(Accessor):
    def __init__(self,path,chrom,sense,system='hg19',offset=0,suffix='.maf.stitched.cmpl.repeats_lc',**kwargs):
        super(MSFAccessor,self).__init__(path,chrom,sense,system=system,**kwargs)

        self.logger = logging.getLogger('byo.io.MSFAccessor')
        self.logger.debug("mmap'ing aligned, stitched genomic sequences for chromosome %s from '%s' offset=%d" % (chrom,path,offset))

        self.system = system
        self.data = None
        self.offset = offset
        
        # try to access the whole genome, using indexing for fast lookup
        trials = [os.path.join(path,chrom+suffix)]
        for fname in trials:
            self.logger.debug("trying to load '%s'" % fname)
            if os.access(fname,os.R_OK):
                self.data = IndexedFasta(fname,**kwargs)
                self.species = sorted(self.data.chrom_stats.keys())
                self.extra_data = dict(species=self.species)
                self.covered_strands = [chrom+'+',chrom+'-']
                
                # mangle the "chrom" names to make sure 
                # reference sequences like 'hg19.chr5(+):1-12345' are found as 'hg19'
                for k,v in self.data.chrom_stats.items():
                    parts = k.split(".chr")
                    if len(parts) > 1:
                        self.data.chrom_stats[parts[0]] = v
                break
                
        if not self.data:
            # all fails: return Ns only
            self.logger.warning("Could not access any of '%s'. Switching to dummy mode (only Ns)" % str(trials))
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy
            self.covered_strands = [chrom+'+',chrom+'-']
            self.species = []
        else:
            # register for all chroms/strands
            self.covered_strands = [chrom+'+',chrom+'-']

        # TODO: maybe remove this if not needed
        self.get = self.get_oriented


    def load_indexed(self,path):
        ipath = path+'.index'
        if not os.access(ipath,os.R_OK):
            index = self.build_index(path,ipath)
        else:
            index = self.load_index(ipath)

        self.chrom_ofs = index.chrom_ofs
        

    def get_data(self,chrom,start,end,sense,species=[]):
        start += self.offset
        end += self.offset

        if not species:
            species = self.species
        if start < 0 or end < 0:
            return [self.get_dummy(spc,start,end,sense) for spc in species]
        
        if sense == "+":
            return [self.data.get_data(spc,start,end,sense).replace('-','N') for spc in species]
        else:
            return [complement(self.data.get_data(spc,start,end,'+').replace('-','N')) for spc in species]


    def get_oriented(self,chrom,start,end,sense,species=[]):
        start += self.offset
        end += self.offset

        if not species:
            species = self.species
        if start < 0 or end < 0:
            return [self.get_dummy(spc,start,end,sense) for spc in species]
        
        if sense == "+":
            return [self.data.get_data(spc,start,end,sense).replace('-','N') for spc in species]
        else:
            return [rev_comp(self.data.get_data(spc,start,end,'+').replace('-','N')) for spc in species]


    def get_dummy(self,chrom,start,end,sense,species=[],**kwargs):
        if not species:
            species = self.species
        
        return ["N"*int(end-start) for spc in species]

        
if __name__ == "__main__":
    #i = IndexedFasta('/data/BIO2/pcp/systems/hg19/genome/hg19.fna')
    i = IndexedFasta('/data/BIO2/pcp/systems/ce6/genome/ce6.fna')
    #print i.chrom_stats
    print i.get_data('chrIV',0,100,'+')
    #print i.chrom_stats
        
