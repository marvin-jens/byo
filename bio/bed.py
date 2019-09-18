# coding=future_fstrings
from byo.bio.lazytables import NamedTupleImporter
import numpy as np
import logging
import time
from collections import defaultdict

bed6_format = "#chrom:str\tstart:int\tend:int\tname:str\tscore:float\tstrand:str"
bed12_format = "#chrom:str\tstart:int\tend:int\tname:str\tscore:float\tstrand:str\tcds_start:int\tcds_end:int\trgb:intlist\tblocks:int\tblock_sizes:intlist\tblock_starts:intlist"

def bed_importer(src, fmt=bed6_format, **kwargs):
    for rec in NamedTupleImporter(src, descr=fmt,parse_comments=False,default_cast='str',keep_line=False,**kwargs):
        yield rec
    

def transcripts_from_bed(src, system=None):
    from byo.gene_model import Transcript
    for bed in src:
        exon_starts = np.array(bed.block_starts) + bed.start
        exon_ends = exon_starts + np.array(bed.block_sizes)
        tx = Transcript(
                bed.name,
                bed.chrom,
                bed.strand,
                exon_starts,
                exon_ends,
                (bed.cds_start, bed.cds_end),
                score=bed.score,
                system=system,
                rgb=bed.rgb
        )
        yield tx


class BEDTrack(object):
    def __init__(self, fname):
        self.fname = fname
        self.logger = logging.getLogger(f'byo.bio.BedTrack({fname})')
        import gzip
        if fname.endswith('.gz'):
            self.bfile = gzip.GzipFile(fname, "r")
            self.logger.debug('Gzip compression detected')
        else:
            self.bfile = file(fname, "r")

        self.starts = defaultdict(list)
        self.ends = defaultdict(list)
        self.ids = defaultdict(list)
        self.bed_records = []
        self.ncls = {}

        self._load_bed()
        self._build_ncls()
    
    def _load_bed(self):
        for line in self.bfile:
            chrom, start, end, name, score, strand = line.rstrip().split('\t')[:6]
            start = int(start)
            end = int(end)
            
            n = len(self.bed_records)
            self.bed_records.append( (chrom, start, end, name, score, strand) )
            strand = chrom+strand

            self.starts[strand].append(start)
            self.ends[strand].append(end)
            self.ids[strand].append(n)

        self.logger.debug('loaded {} BED records'.format(len(self.bed_records)))

    def _build_ncls(self):
        import ncls
        t0 = time.time()
        for strand in sorted(self.ids.keys()):
            starts = np.array(self.starts[strand])
            ends = np.array(self.ends[strand])
            ids = np.array(self.ids[strand])
            # print strand, starts, ends, ids
            self.ncls[strand] = ncls.NCLS(
                starts,
                ends,
                ids,
            )

        self.bed_records = np.array(self.bed_records, dtype=object)
        dt = time.time() - t0
        self.logger.debug("built NCLists for {0} items in {1:.3f}seconds".format(len(self.bed_records), dt))
    
    def get(self, chrom, start, end, sense):
        strand = chrom+sense
        ncls = self.ncls[strand]
        found = [self.bed_records[fid[2]] for fid in self.ncls[strand].find_overlap(start, end)]

        return found

    def get_oriented(self, chrom, start, end, sense):
        return self.get(chrom, start, end, sense)

