import sys,re
import numpy
from bisect import bisect_left

def MAF_chunks(source,ref,species=[]):

    def var_parse(line):
        d = {}
        for m in re.finditer(r'(\S+)=(\S+)',line):
            k,v = m.groups()
            d[k] = v
        return d
        
    maf_dict = {}
    comments = []

    class MAFBlock(object):
        def __init__(self):
            self.maf_info = maf_dict
            self.comments = comments
            self.seqs = {}
            self.info = {}
            self.species = []
            self.spec_coords = {}
            self.chrom = "??"
            self.start = -1
            self.end = -1
            self.size = 0
            self.lines = []
            self.n_cols = 0
            #print "init",len(self.species)
        def _setup(self):
            #print "setup",len(self.species)
            if self.size > 0:
                return

            seq = self.seqs[ref]
                
            block_starts = []
            block_ends = []
            
            gap_start = -1
            if seq[0] != '-':
                block_starts.append(0)
            else:
                gap_start = 0

            for i,l in enumerate(seq):
                if gap_start >=0 and l != "-":
                    block_starts.append(i)
                    gap_start = -1
                elif gap_start == -1 and l == "-":
                    gap_start = i
                    block_ends.append(i)
                if l != '-':
                    self.size += 1

            if seq[-1] != '-':
                block_ends.append(len(seq))

            self.block_starts = numpy.array(block_starts)
            self.block_ends = numpy.array(block_ends)
            self.block_lengths = self.block_ends - self.block_starts
            self.ref_starts = numpy.array([0,]+ list(self.block_lengths.cumsum()))
            
            #print zip(self.block_starts,self.block_ends)
            #print self.ref_starts
        
        @property
        def maf_str(self):
            return "\n".join(self.lines)
            
        def add_sequence(self,species,chrom,start,end,strand,seq):
            #print "add_sequence",species,len(self.species)
            if species == ref:
                self.chrom = chrom
                self.start = start
                self.end = end                
                self.strand = strand

            self.seqs[species] = seq
            self.species.append(species)
            self.spec_coords[species] = (chrom,start,end,strand)
            setattr(self,species,seq)

        def get_windows(self,pos,flank=20):
            mpos = self.ref_to_maf(pos)
            windows = {}

            for s,seq in self.seqs.items():
                
                # walk right and collect nucleotides
                win = []
                x = mpos
                while len(win) < flank+1:
                    if x >= len(seq):
                        win.append('N')

                    elif seq[x] != '-':
                        win.append(seq[x])

                    x += 1

                # walk left
                x = mpos
                while len(win) < 2*flank+1:
                    if x < 0:
                        win.insert(0,'N')
                    
                    elif seq[x] != '-':
                        win.insert(0,seq[x])
                        
                    x -= 1

                #print "%20s" % s,seq[mpos-10:mpos+10]
                windows[s] = "".join(win)
                #print "%20s" % s,"".join(win)

            return windows
        
        @property
        def columns(self):
            species = self.species #sorted(self.seqs.keys())
            matrix = [list(self.seqs[s]) for s in species]
            matrix = numpy.array(matrix).transpose()

            return matrix,species
            
        def ref_to_maf(self,pos):
            self._setup()

            p = pos - self.start
            #print pos,self.start,p
            
            n = bisect_left(self.ref_starts,p+1) -1
            #print n
            return self.block_starts[n] + p - self.ref_starts[n]
            
        def __str__(self):
            block = ["MAFBlock(%s%s:%d-%d)" % (self.chrom,self.strand,self.start,self.end),str(self.maf_info)] + self.comments
            for s in self.species:
                block.append("%20s %s" % (s,self.seqs.get(s,'='*self.n_cols)))

            return "\n".join(block)

        def __len__(self):
            return self.n_cols

    block = MAFBlock()#(species=species)

    for line in source:
        #print "MAFLINE",line
        line = line.rstrip()
        block.lines.append(line)
        
        if line.startswith("##"):
            maf_dict.update(var_parse(line))

        elif line.startswith('#'):
            comments.append(line)

        elif line.startswith('a'):
            # alignment block description line
            if len(block):
                # yield preceding block!
                yield block
                block = MAFBlock()
            # keep the description of the next one
            block.block_info = var_parse(line)
            
        elif line.startswith('s'):
            # record the sequence information
            parts = re.split(r'\s+',line.strip())
            
            species,chrom = re.match(r'(\w+)\.(\S+)',parts[1]).groups()
            start,size,strand,srcSize = parts[2:6]
            
            if strand == "+":
                start = int(start)
                end = start + int(size)
                
            else:
                end = int(srcSize) - int(start)
                start = end - int(size)
            
            block.add_sequence(species,chrom,start,end,strand,parts[-1])
            block.n_cols = len(parts[-1])
                
        elif line.startswith('i'):
            # record the information relating to neighboring blocks
            parts = re.split(r'\s+',line.strip())
            #print parts
            species,chrom = re.match(r'(\w+)\.(\S+)',parts[1]).groups()
            if len(parts[2:]) != 4:
                sys.stderr.write('malformed line skipped:"%s"\n' % line.rstrip())
                continue

            prv_code,prv_num,nxt_code,nxt_num = parts[2:]
            block.info[species] = (prv_code,int(prv_num),nxt_code,int(nxt_num))
        else:
            # line is empty or broken or whatever
            if line.strip():
                sys.stderr.write('malformed, unhandled line skipped:"%s"\n' % line.rstrip())
            continue
        
    if len(block):
        yield block
        block = MAFBlock()

if __name__ == "__main__":
    from gzip import GzipFile
    for i,maf in enumerate(MAF_chunks(GzipFile('/data/rajewsky/genomes/hg19_100way/maf/chr22.maf.gz'),'hg19')):
        print maf
        maf._setup()
        center = maf.start+maf.size/2
        #cmaf = maf.ref_to_maf(center)
        #print maf.size,center,maf.size/2,cmaf

        print "windows"
        print maf.get_windows(center)
        
        print "columns"
        print maf.columns
        if i > 10:
            break