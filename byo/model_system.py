import os,re, sys

from byo.track import Track, load_track
from byo.bio.genome_accessor import GenomeCache, RemoteCache
from byo.bio.annotation import AnnotationAccessor
import byo.config
import logging

class LazyTranscriptLoader(object):
    def __init__(self,system = None):
        self.transcripts = None
        self.system = system
        self.logger = logging.getLogger("LazyTranscriptLoader(system={system.name})".format(system=self.system) )

    def __getitem__(self,txname):
        if not self.transcripts:
            self.transcripts = self.load_transcript_catalogs()

        return self.transcripts[txname]

    def load_transcript_catalogs(self):
        from byo.gene_model import transcripts_from_UCSC
        import glob
        path = os.path.join(self.system.root, "annotation", self.system.name, "*.ucsc.gz")
        sources = glob.glob(path)
        if sources:
            self.logger.debug('loading {0}'.format(sources[0]))
            T = transcripts_from_UCSC(sources[0],system = self.system)
            for s in sources[1:]:
                self.logger.debug('loading {0}'.format(s))
                T.load(s)

            self.logger.info("loaded {0} transcript models from {1} source(s)".format(len(T), len(sources)))
            return T
        else:
            self.logger.error("no transcript models found in path '{0}'".format(path))
            return {}

class ModelSystem(object):
    def __init__(self, name, genome = None, transcript_models = None, root = byo.config.system_root):
        self.name = name
        self.root = root
        if genome == None:
            if getattr(byo.config, "genome_server", None):
                # print "getting remote genome from", byo.config.genome_server, "for", name
                self.genome = RemoteCache(byo.config.genome_server)[name]
            else:
                self.genome = GenomeCache(os.path.join(root, "genomes"))[name]
        else:
            self.genome = genome

        if transcript_models == None:
            self.transcript_models = LazyTranscriptLoader(system=self)
        else:
            self.transcript_models = transcript_models
            
        # Fails before first access due to lazy loading of genome
        #self.chr_sizes = self.genome.data.chrom_stats
        

    def get_annotations_track(path="",accessor=AnnotationAccessor,**kwargs):
        if not path:
            path = os.path.join(self.root,"annotation",self.name,"compiled")
        return Track(path,accessor,**kwargs)

    def load_track(self, path, **kwargs):
        return load_track(path, system=self, **kwargs)

    def get_refGenes(self):
        return self.transcript_models
