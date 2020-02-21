# -*- coding: future_fstrings -*-
from __future__ import print_function
import sys
import os
try:
    import cPickle as pickle
except ImportError:
    import pickle
    
import time
from collections import defaultdict
import byo
import byo.gene_model
import byo.annotation
import numpy as np
import pandas as pd

import logging
logging.basicConfig(level=logging.DEBUG)
    
ann = byo.annotation.AnnotationTrack('hg19/gencode28', url="tcp://heart-of-gold.mit.edu:13370")

df = pd.read_table(sys.stdin, header=None, names=["chrom", "start", "end", "name", "score", "sense", ])
df.sort_values("score", ascending=False, inplace=True)

genome = byo.systems['hg19'].genome

for x in df.itertuples():
    q = ann.get_oriented(x.chrom, x.start, x.end, x.sense)
    seq = genome.get_oriented(x.chrom, x.start-20, x.end+20, x.sense)

    cat_counts, tx_set, tx_exon_set = byo.annotation.summarize(q)
    # print(q.match_identifiers)
    # print(cat_counts)
    kind = byo.annotation.categorize(cat_counts)

    out = [x.chrom, x.start, x.end, x.name, x.score, x.sense, kind, ",".join(sorted(tx_set)), ",".join(sorted(tx_exon_set)), seq]
    print("\t".join([str(o) for o in out]))

