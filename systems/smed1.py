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

genome = Track(os.path.join(root,"smed1","genome"),GenomeAccessor,system='smed1')

def get_annotation_track(path=os.path.join(root,"smed1","system_annotation"),**kwargs):
    return Track(path,AnnotationAccessor,**kwargs)
    
chr_sizes = {}
try:
    for c in Importer(os.path.join(root,"smed1","chrom.sizes"),descr="## chrom:str \t size:int"):
        chr_sizes[c.chrom] = c.size
except IOError:
    loaded = False
else:
    loaded = True

    

def get_refGenes(path = [os.path.join(root,'smed1','annotation','dresden.ucsc'), os.path.join(root,'smed1','annotation','bimsb.ucsc')]):
    from sequence_data.gene_model import transcripts_from_UCSC
    gene_names = {}
    #for l in file(os.path.join(root,'annotation','gene_names')):
    #    k,v = l.split('\t')
    #    gene_names[k] = v.rstrip()

    import sequence_data
    #return transcripts_from_UCSC(path,system = sequence_data.smed1,gene_names=gene_names,fix_chr=False)
    if type(path) == str:
        return transcripts_from_UCSC(path,system=sequence_data.smed1,fix_chr=False)
    else:
        T = transcripts_from_UCSC(path[0],system=sequence_data.smed1,fix_chr=False)
        for p in path[1:]:
            T.load(p)
        return T

