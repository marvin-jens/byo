# -*- coding: utf-8 -*-
from byo.track import Track
from byo.io.track_accessors import Accessor, ArrayAccessor, GenomeAccessor, AnnotationAccessor
from byo.io.lazytables import NamedTupleImporter as Importer
import byo.config

import os,re
from logging import debug,warning,error,info
from numpy import uint32, float32

root = byo.config.system_root

genome = Track(os.path.join(root,"reference","hg19"),GenomeAccessor,system='hg19')

def get_annotation_track(path=os.path.join(root,"annotation","hg19","compiled"),accessor=AnnotationAccessor,**kwargs):
    return Track(path,accessor,**kwargs)

def get_refGenes(path = [
    os.path.join(root,'annotation','hg19','wgEncodeGencodeBasicV17.ucsc'),
    #os.path.join(root,'annotation','hg19','refGene.ucsc'),
    #os.path.join(root,'annotation','hg19','transcript.noncoding.lincRNA.ucsc_lincrna_track.ucsc')
    ] ):

    from byo.gene_model import transcripts_from_UCSC
    import sys
    if type(path) == str:
        return transcripts_from_UCSC(path,system = sys.modules[__name__])
    else:
        T = transcripts_from_UCSC(path[0],system = sys.modules[__name__])
        for p in path[1:]:
            T.load(p)
        return T

chr_sizes = {}
try:
    for c in Importer(os.path.join(root,'reference','hg19',"chrom.sizes"),descr="## chrom:str \t size:int"):
        chr_sizes[c.chrom] = c.size
except IOError:
    loaded = False
else:
    loaded = True
