
import numpy as np
import bx.bbi
from bx.bbi.bigwig_file import BigWigFile
import logging

class BigWigTrack(object):
    def __init__(self, fplus, fminus=None, minus_flip_sign=False):
        self.logger = logging.getLogger('byo.bio.BigWigTrack(+={fplus} -={fminus})')
        plus = BigWigFile(file(fplus, 'rb'))
        if plus is None:
            self.logger.error("BW plus not found")
        
        minus = BigWigFile(file(fminus, 'rb'))
        if minus is None:
            self.logger.error("BW minus not found")
        
        if not (plus is None) and not (minus is None):
            self.stranded = True
            self.strands = {
                '+' : plus,
                '-' : minus,
            }
        else:
            self.stranded = False
            self.strands = {
                '+' : plus,
                '-' : plus,
            }
        
        self.minus_flip_sign = minus_flip_sign

    def get(self, chrom, start, end, strand):
        bw = self.strands[strand]
        data = bw.get_as_array(chrom, start, end)
        if strand == '-' and self.minus_flip_sign:
            data = - data # undo the minus sign included in e.g. ENCODE eCLIP tracks
        # print f"reading data for {chrom}:{start}-{end}{strand} -> {data}"
        return data

    def get_oriented(self, chrom, start, end, strand):
        data = self.get(chrom, start, end, strand)
        if strand == '-':
            data = data[::-1]
        return data
