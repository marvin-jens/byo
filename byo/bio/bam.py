# coding=future_fstrings
import pysam
import os
import itertools
import numpy as np
import logging

class BAMTrack(object):
    def __init__(self, fname, mode='list', only_once=False):
        self.fname = fname
        self.bam = pysam.AlignmentFile(self.fname)
        self.mode = mode
        self.logger = logging.getLogger(f'byo.BAMTrack({fname})')
        try:
            has_index = self.bam.check_index()
        except ValueError:
            self.logger.debug("check_index() raised exception")
            has_index = False
    
        if not has_index:
            self.logger.info(f"indexing BAM file... {fname}")
            pysam.samtools.index(fname)
            # re-open BAM file, otherwise it still thinks index is missing :/
            self.bam = pysam.AlignmentFile(self.fname)

        get_functions = {
            'list' : (self._get_list, self._get_list_oriented),
            '5p_array' : (self._get_5p, self._get_oriented),
            '3p_array' : (self._get_3p, self._get_oriented),
            'coverage' : (self._get_cov, self._get_oriented),
            'count' : (self._get_count, self._get_count),
        }

        self.get, self.get_oriented = get_functions[mode]
        self.logger.debug(f"initialized track in '{mode}' mode")
        self.seen = set()
        self.only_once = only_once

    @property
    def total_mapped_reads(self):
        return self.bam.mapped

    def _get(self, chrom, start, end, strand, mate=None):
        def strand_filter(read):
            if strand == '*':
                return True
            
            rev = read.is_reverse
            if strand == '+' and not rev:
                return True
                
            elif strand == '-' and rev:
                return True

            else:
                return False

        def mate_filter(read):
            if mate is None:
                return True
            
            m1 = read.is_read1
            m2 = read.is_read2
            if mate == 1 and m1:
                return True
            
            elif mate == 2 and m2:
                return True
            
            else:
                return False

        def once_filter(read):
            if read.query_name in self.seen:
                return False

            self.seen.add(read.query_name)
            return True

        reads = []
        for read in self.bam.fetch(chrom, start, end):
            if not strand_filter(read):
                continue

            if self.only_once and not once_filter(read):
                continue

            if not mate_filter(read):
                continue

            reads.append(read)
        
        return reads

    def _get_count(self, *argc, **kw):
        reads = self._get(*argc, **kw)
        return len(reads)

    def _get_list(self, chrom, start, end, strand, **kw):
        reads = self._get(chrom, start, end, strand, **kw)
        return sorted(reads, key=lambda r : r.pos)

    def _get_list_oriented(self, chrom, start, end, strand, **kw):
        reads = self._get(chrom, start, end, strand, **kw)
        if strand == '-':
            return sorted(reads, key=lambda r : r.aend)
        else:
            return sorted(reads, key=lambda r : r.pos)

    def _get_5p(self, chrom, start, end, strand, **kw):
        reads = self._get(chrom, start, end, strand, **kw)
        res = np.zeros(end-start)
        for r in reads:
            x = r.pos if strand == '+' else r.aend
            if start <= x < end:
                res[x-start] += 1
        
        return res

    def _get_3p(self, chrom, start, end, strand, **kw):
        reads = self._get(chrom, start, end, strand, **kw)
        res = np.zeros(end-start)
        for r in reads:
            x = r.pos if strand == '-' else r.aend
            if start <= x < end:
                res[x-start] += 1
        
        return res

    def _get_cov(self, chrom, start, end, strand, **kw):
        reads = self._get(chrom, start, end, strand, **kw)
        res = np.zeros(end-start)
        for r in reads:
            x, y = r.pos, r.aend
            s = max(0, x - start)
            e = min(end - start, y - start)
            res[s:e] += 1

        return res

    def _get_oriented(self, chrom, start, end, strand, **kw):
        res = self.get(chrom, start, end, strand, **kw)
        if strand == '-':
            res = res[::-1]
        
        return res