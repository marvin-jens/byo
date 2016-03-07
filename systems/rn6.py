# -*- coding: utf-8 -*-
from byo.track import Track
from byo.io.track_accessors import Accessor, ArrayAccessor, GenomeAccessor ,AnnotationAccessor
from byo.io.lazytables import NamedTupleImporter as Importer
from byo.io.lazytables import LazyImporter
import byo.config
import byo
import os,re

from logging import debug,warning,error,info

root =  byo.config.system_root

genome = Track(os.path.join(root,"rn6","genome"),GenomeAccessor,system='rn6')

def get_annotation_track(path=os.path.join(root,"rn6","system_annotation"),**kwargs):
    return Track(path,AnnotationAccessor,**kwargs)

#def get_refGenes(path = os.path.join(root,'annotation','flyBaseGene.ucsc')):
def get_refGenes(path = [
    os.path.join(root,"rn6",'annotation','ensGene.ucsc'),
    os.path.join(root,"rn6",'annotation','sgpGene.ucsc'),
    ] ):
    from byo.gene_model import transcripts_from_UCSC
    import sys
    
    gene_names = {}
    for l in file(os.path.join(root,"rn6",'annotation','gene_names')):
        k,v = l.split('\t')
        gene_names[k] = v.rstrip()

    if type(path) == str:
        return transcripts_from_UCSC(path,system = sys.modules[__name__],gene_names=gene_names)
    else:
        T = transcripts_from_UCSC(path[0],system = sys.modules[__name__],gene_names=gene_names)
        for p in path[1:]:
            T.load(p)
        return T

    
chr_sizes = {}
try:
    for c in Importer(os.path.join(root,"rn6","chrom.sizes"),skip=[2],descr="## chrom:str \t size:int \t fileName:str"):
        chr_sizes[c.chrom] = c.size
except IOError:
    loaded = False
else:
    loaded = True


    