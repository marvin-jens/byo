# coding=future_fstrings
from byo.bio.lazytables import NamedTupleImporter
import numpy as np
import logging
import time
from collections import defaultdict, namedtuple

bed6_format = "#chrom:str\tstart:int\tend:int\tname:str\tscore:float\tstrand:str"
bed12_format = "#chrom:str\tstart:int\tend:int\tname:str\tscore:float\tstrand:str\tcds_start:int\tcds_end:int\trgb:intlist\tblocks:int\tblock_sizes:intlist\tblock_starts:intlist"

bed6_tuple = namedtuple("BED6", 'chrom start end name score strand')
bed12_tuple = namedtuple("BED12", 'chrom start end name score strand cds_start cds_end rgb blocks block_sizes block_starts')
def bed_importer(src, fmt=bed6_format, **kwargs):
    for rec in NamedTupleImporter(src, descr=fmt,parse_comments=False,default_cast='str',keep_line=False,**kwargs):
        yield rec
    

def transcript_from_bed6(bed, system=None):
    exon_starts = np.array([bed.start, ])
    exon_ends = np.array([bed.end, ])
    tx = Transcript(
            bed.name,
            bed.chrom,
            bed.strand,
            exon_starts,
            exon_ends,
            (bed.start, bed.start),
            score=bed.score,
            system=system,
    )
    return tx

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
    def __init__(self, fname, bed_type=None, system=None, keep_extra_columns=False, row_filter=None, pad_us=0, pad_ds=0, **kw):
        self.fname = fname
        self.system = system
        self.logger = logging.getLogger(f'byo.bio.BedTrack({fname})')
        import gzip
        if fname.endswith('.gz'):
            self.bfile = gzip.GzipFile(fname, "r")
            self.logger.debug('Gzip compression detected')
        else:
            self.bfile = file(fname, "r")

        self.pad_us = pad_us
        self.pad_ds = pad_ds
        self.logger.debug(f'pad_us={self.pad_us} pad_ds={self.pad_ds}')
        self.starts = defaultdict(list)
        self.ends = defaultdict(list)
        self.ids = defaultdict(list)
        self.bed_records = []
        self.ncls = {}

        self._bed_type = bed_type
        self.keep_extra_columns = keep_extra_columns
        self.n_cols = 6
        self.row_filter = row_filter
        self._load_bed()
        self._build_ncls()
    
    def _load_bed(self):
        bed_tuple = bed6_tuple
        for line in self.bfile:
            if self._bed_type is None:
                # determine BED format from column count:
                pass # TODO: needs implementation
                
            cols = line.rstrip().split('\t')
            if not self.row_filter is None:
                if not self.row_filter(cols):
                    continue

            chrom, start, end, name, score, strand = cols[:6]
            start = int(start)
            end = int(end)
            if strand == '+':
                start -= self.pad_us
                end += self.pad_ds
            else:
                start -= self.pad_ds
                end += self.pad_us

            assert(end - start >= (self.pad_us + self.pad_ds))
            cols[1:3] = [start, end]
            if self.keep_extra_columns and len(cols) > 6 and bed_tuple == bed6_tuple:
                fmt = 'chrom start end name score strand ' + \
                    " ".join([f'extra{i}' for i in range(len(cols)-6)])
                # print(fmt)
                bed_tuple = namedtuple('BED6_plus', fmt)
                self.n_cols = len(cols)

            cols = cols[:self.n_cols]

            n = len(self.bed_records)
            # print(len(cols), cols)
            self.bed_records.append( bed_tuple(*cols) )
            strand = chrom+strand

            self.starts[strand].append(start)
            self.ends[strand].append(end)
            self.ids[strand].append(n)

        n = len(self.bed_records)
        self.logger.debug(f'loaded {n} BED records')

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

        self.bed_records = np.array(self.bed_records + [None], dtype=object)
        # None is added bc it prevents the namedtuples to convert to pure numpy


        dt = time.time() - t0
        self.logger.debug("built NCLists for {0} items in {1:.3f}seconds".format(len(self.bed_records), dt))
    
    def get(self, chrom, start, end, sense):
        strand = chrom+sense
        if strand in self.ncls:
            ncls = self.ncls[strand]
            found = [self.bed_records[fid[2]] for fid in self.ncls[strand].find_overlap(start, end)]
        else:
            found = []
        return found

    def get_oriented(self, chrom, start, end, sense):
        return self.get(chrom, start, end, sense)

    def __iter__(self):
        from byo.gene_model import Transcript

        for bed in self.bed_records:
            if bed is None:
                return
            tx = Transcript(bed.name, bed.chrom, bed.strand, [bed.start,], [bed.end,], (bed.start, bed.start), score=bed.score, system=self.system)
            tx._bed = bed
            yield tx

