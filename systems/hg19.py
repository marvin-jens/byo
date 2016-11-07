# -*- coding: utf-8 -*-
from byo.track import Track
from byo.io.track_accessors import Accessor, ArrayAccessor, GenomeAccessor, AnnotationAccessor
from byo.io.lazytables import NamedTupleImporter as Importer
import byo.config

import os,re, sys
from logging import debug,warning,error,info
from numpy import uint32, float32
from byo.io.genome_accessor import GenomeCache

root = byo.config.system_root
name = "hg19"
genome = GenomeCache(os.path.join(root,"genomes"))[name]

class LazyTranscriptLoader(object):
    def __init__(self,system = None):
        self.transcripts = None
        self.system = system

    def __getitem__(self,txname):
        if not self.transcripts:
            self.transcripts = self.system.get_refGenes()
        return self.transcripts[txname]

transcript_models = LazyTranscriptLoader(system = sys.modules[__name__])

def get_annotation_track(path=os.path.join(root,"annotation",name,"compiled"),accessor=AnnotationAccessor,**kwargs):
    return Track(path,accessor,**kwargs)

def get_refGenes(path = [
    os.path.join(root,'annotation',name,'wgEncodeGencodeBasicV17.ucsc'),
    #os.path.join(root,'annotation',name,'refGene.ucsc'),
    #os.path.join(root,'annotation',name,'transcript.noncoding.lincRNA.ucsc_lincrna_track.ucsc')
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
    for c in Importer(os.path.join(root,'reference',name,"chrom.sizes"),descr="## chrom:str \t size:int"):
        chr_sizes[c.chrom] = c.size
except IOError:
    loaded = False
else:
    loaded = True
