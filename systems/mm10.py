# -*- coding: utf-8 -*-
from byo.track import Track
from byo.io.track_accessors import Accessor, ArrayAccessor, GenomeAccessor ,AnnotationAccessor
from byo.io.lazytables import NamedTupleImporter as Importer
from byo.io.lazytables import LazyImporter
import byo.config
import byo
import os,re, sys

from logging import debug,warning,error,info
from byo.io.genome_accessor import GenomeCache

name = "mm10"

root =  byo.config.system_root
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


#def get_refGenes(path = os.path.join(root,'annotation','flyBaseGene.ucsc')):
def get_refGenes(path = os.path.join(root,"mm10",'annotation','ensGene.ucsc')):
    from byo.gene_model import transcripts_from_UCSC
    gene_names = {}
    for l in file(os.path.join(root,"mm10",'annotation','gene_names')):
        k,v = l.split('\t')
        gene_names[k] = v.rstrip()

    import sys
    return transcripts_from_UCSC(path,system = sys.modules[__name__],gene_names=gene_names)
    
chr_sizes = {}
try:
    for c in Importer(os.path.join(root,'reference',name,"chrom.sizes"),descr="## chrom:str \t size:int"):
        chr_sizes[c.chrom] = c.size
except IOError:
    loaded = False
else:
    loaded = True
   
