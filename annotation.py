import sys
import time
import re
import logging
import ncls
import zmq
import byo
import byo.gene_model
import numpy as np
from collections import defaultdict

terminals = set(['PAS', "5'SS", "3'SS", "CDS", "5UTR", "3UTR"])

def categorize(cats):
    res = ''
    n = float(np.array(cats.values()).max())
    exon_frac = cats['exon'] / n
    intron_frac = cats['intron'] / n
    if exon_frac > .1:
        res = 'exon'
    else:
        return 'intron'
        # res = 'intron'
    
    segs = ["5UTR", "CDS", "3UTR"]
    fracs = np.array([cats[s]/ n for s in segs])

    seg = segs[fracs.argmax()]

    return seg + "_" + res

def cat_from_name(name):
    parts = name.split('/')
    end = parts[-1]
    if end not in terminals:
        if end.startswith('exon'):
            return parts[0], 'exon'
        
        elif end.startswith('intron'):
            return parts[0], 'intron'

    return parts[0], end

def summarize(q):
    # print q
    cat_count = defaultdict(int)
    tx_set = set()
    tx_exon_set = set()
    for name in q.match_identifiers:
        tx_id, cat = cat_from_name(name)
        tx_set.add(tx_id)
        cat_count[cat] += 1
        if cat == 'exon':
            tx_exon_set.add(tx_id)
    
    return cat_count, tx_set, tx_exon_set


class Query(object):
    def __init__(self, realm, coord=None, collapsed=False, unique=False, identifier=None, get_objects=False, **kwargs):
        self.realm = realm
        self.coord = coord
        self.collapsed = collapsed
        self.unique = unique
        if unique:
            # unique implies collapsed or things get inconsistent
            self.collapsed = True

        self.get_objects = get_objects
        self.identifier = identifier

    def __str__(self):
        s = "Query({self.realm}, coord={self.coord}, identifier={self.identifier}, collapsed={self.collapsed}, unique={self.unique})".format(self=self)
        if hasattr(self, "status"):
            s += " status={self.status}".format(self=self)
        if hasattr(self, "match_identifiers"):
            s += "-> matches={}".format(",".join(self.match_identifiers))
        
        return s
    
    def __add__(self, q):
        if self.realm != q.realm:
            raise ValueError("realm mismatch!")
        
        Q = Query(self.realm)

        # concatenate matches
        Q.match_starts = getattr(self, "match_starts", []) + getattr(q, "match_starts", [])
        Q.match_ends = getattr(self, "match_ends", []) + getattr(q, "match_ends", [])
        Q.match_identifiers = getattr(self, "match_identifiers", []) + getattr(q, "match_identifiers", [])
        Q.match_objects = getattr(self, "match_objects", []) + getattr(q, "match_objects", [])

        return Q
 

class AnnotationTrack(object):
    def __init__(self, realm, url="tcp://*:13370"):
        self.realm = realm
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.REQ)
        self.url = url
        self.socket.connect(self.url)

    def get_oriented(self, chrom, start, end, sense, **kwargs):
        query = Query(self.realm, coord = (chrom, int(start), int(end), sense), **kwargs)
        self.socket.send_pyobj(query)
        res = self.socket.recv_pyobj()
        return res

    def get_transcripts_by_id(self, transcript_id, **kwargs):
        query = Query(self.realm, identifier=("transcript_id", transcript_id), **kwargs)
        self.socket.send_pyobj(query)
        res = self.socket.recv_pyobj()
        return res

    def get_transcripts_by_gene(self, gene_id, **kwargs):
        query = Query(self.realm, identifier=("gene_id", gene_id), **kwargs)
        self.socket.send_pyobj(query)
        res = self.socket.recv_pyobj()
        return res


class ExpressionTrack(object):
    def __init__(self, fexpr, ann_track, pattern="exon"):
        self.ann_track = ann_track
        self.pattern = pattern
        self.expr, self.expr_err = load_expression(fexpr)
    
    def combined_expression(self, chain):
        def list_join(ll):
            x = []
            for l in ll:
                x.extend(l)
            return x

        tpm = 0
        for start, end, (feature_id, tx_model) in chain.splice(self.ann_track, join=list_join, pattern=self.pattern):
            name = feature_id.split('.')[0]
            # print name, coll.expr[name]
            tpm += coll.expr[name]
        return tpm

           
class AnnotationServer(object):
    def __init__(self, realm="annotation"):
        self.realm = realm
        self.starts = defaultdict(list)
        self.ends = defaultdict(list)
        self.ids = defaultdict(list)
        self.features = defaultdict(list)
        self.n_feat = 0
        self.n_tx = 0
        self.ncls = {}
        self.logger = logging.getLogger("byo.AnnotationServer")
        self.transcripts = {}

        self.tx_by_id = defaultdict(list)
        self.tx_by_gene = defaultdict(list)

    def build_ncls(self):
        for strand in sorted(self.ids.keys()):
            t0 = time.time()
            starts = np.array(self.starts[strand])
            ends = np.array(self.ends[strand])
            ids = np.array(self.ids[strand])
            # print strand, starts, ends, ids
            self.ncls[strand] = ncls.NCLS(
                starts,
                ends,
                ids,
            )
            self.features[strand] = np.array(self.features[strand])
            dt = time.time() - t0
            self.logger.debug("built NCList from {0} items in {1:.3f}seconds".format(len(starts), dt))

    def add_transcript(self, tx, feature_extractor = lambda tx : list(tx.features) + list(tx.segments)):
        self.transcripts[tx.transcript_id] = tx
        self.tx_by_id[tx.name.split('.')[0]].append(tx)
        self.tx_by_gene[tx.gene_id].append(tx)
        self.tx_by_gene[tx.gene_name].append(tx)
        self.n_tx += 1
        strand = tx.chrom + tx.sense
        to_add = feature_extractor(tx)
        # print "feat extr", feature_extractor
        for f in to_add:
            # print f
            self.starts[strand].append(f.start)
            self.ends[strand].append(f.end)
            self.ids[strand].append(len(self.features[strand]))
            self.features[strand].append(f.name)
            self.n_feat += 1

    def load_transcripts(self, fname, system='hg19', **kwargs):
        t0 = time.time()
        for tx in byo.gene_model.transcripts_from_GTF(fname, system=system):
            self.add_transcript(tx, **kwargs)
            if not self.n_tx % 1000:
                dt = time.time() - t0
                self.logger.debug("{0} transcripts loaded, {1:.1f} tx/sec".format(self.n_tx, self.n_tx / dt))

        dt = time.time() - t0
        self.logger.debug("loaded {0} transcripts with {1} features in {2:.3f}seconds".format(self.n_tx, self.n_feat, dt))

    def query(self, chrom, start, end, sense):
        strand = chrom+sense
        features = [(fid[0], fid[1], self.features[strand][fid[2]]) for fid in self.ncls[strand].find_overlap(start, end)]
        return features

    def Q(self, q):
        if self.realm != q.realm:
            q.status = (1, "realm mismatch {} != {}".format(self.realm, q.realm))
            return q
        
        if q.coord:
            # coordinate based query!
            chrom, start, end, sense = q.coord
            strand = chrom+sense
            nc = self.ncls[strand]
            fs = self.features[strand]
            if q.collapsed:
                # we only care about *which* features we overlap, not where exactly. 
                # This is much faster, especially for long ranges!
                bla, fids = nc.all_overlaps_both(
                    np.array([start,]), 
                    np.array([end,]), 
                    np.zeros(1, dtype=int)
                )
            else:
                # this reports the overlapping feature start and end coordinates as well
                q.match_starts, q.match_ends, fids = np.array(list(nc.find_overlap(start, end))).T
            if q.unique:
                fids = sorted(set(fids))

            q.match_identifiers = fs[fids]
            if q.get_objects:
                q.match_objects = [self.transcripts[fid.split('/')[0]] for fid in q.match_identifiers]
        
        elif q.identifier:
            kind, name = q.identifier
            if kind == "gene_id" and name in self.tx_by_gene:
                q.match_objects = self.tx_by_gene[name]

            elif kind == "transcript_id" and name in self.tx_by_id:
                q.match_objects = self.tx_by_id[name]

            else:
                q.status = (1, "unknown {} '{}'".format(kind, name))
                return q

            q.match_identifiers = [tx.name for tx in q.match_objects]

        q.status = (0, "OK")
        return q

    def serve_forever(self, bind_addr="tcp://*:13370"):
        context = zmq.Context()
        socket = context.socket(zmq.REP)
        socket.bind(bind_addr)

        while True:
            query = socket.recv_pyobj()
            # print "received request", query
            res = self.Q(query)
            # print "sending result", res
            socket.send_pyobj(res)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    ann = AnnotationServer(realm=sys.argv[1])
    ann.load_transcripts(sys.argv[2])
    ann.build_ncls()

    # q = ann.Q(Query("hg19/gencode28", identifier=("gene_id", "SAMD11")))
    # exon = list(q.match_objects[0].exons)[1]
    # print exon
    # q.identifier = None
    # q.coord = (exon.chrom, exon.start, exon.end, exon.sense)
    # q.collapsed = True
    # q.unique = True
    # q = ann.Q(q)
    # print q
    # for start, end, name in zip(q.match_starts, q.match_ends, q.match_identifiers):
    #     print start, end, name

    ann.serve_forever(bind_addr="tcp://*:13370")

