# -*- coding: utf-8 -*-
from byo.track import Track
from byo.io.track_accessors import Accessor, ArrayAccessor, GenomeAccessor ,AnnotationAccessor
from byo.io.lazytables import NamedTupleImporter as Importer
from byo.io.lazytables import LazyImporter
import byo.config
import byo
import os,re,sys

from logging import debug,warning,error,info

name = "dm3"
root =  byo.config.system_root

genome = Track(os.path.join(root,name,"genome"),GenomeAccessor,system='dm3')

def get_annotation_track(path=os.path.join(root,name,"system_annotation"),**kwargs):
    return Track(path,AnnotationAccessor,**kwargs)


class LazyTranscriptLoader(object):
    def __init__(self,system = None):
        self.transcripts = None
        self.system = system

    def __getitem__(self,txname):
        if not self.transcripts:
            self.transcripts = self.system.get_refGenes()
        return self.transcripts[txname]

transcript_models = LazyTranscriptLoader(system = sys.modules[__name__])

#def get_refGenes(path = os.path.join(root,'annotation','flyBaseGene.ucsc')):
def get_refGenes(path = os.path.join(root,name,'annotation','ensGene.ucsc')):
    from byo.gene_model import transcripts_from_UCSC
    gene_names = {}
    for l in file(os.path.join(root,name,'annotation','gene_names')):
        k,v = l.split('\t')
        gene_names[k] = v.rstrip()

    import sys
    return transcripts_from_UCSC(path,system = sys.modules[__name__],gene_names=gene_names)
    
chr_sizes = {}
try:
    for c in Importer(os.path.join(root,name,"chrom.sizes"),skip=[2],descr="## chrom:str \t size:int \t fileName:str"):
        chr_sizes[c.chrom] = c.size
except IOError:
    loaded = False
else:
    loaded = True


    