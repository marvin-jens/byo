from byo.io.lazytables import NamedTupleImporter

def bed_importer(src,**kwargs):
    bed6_default = "#chrom:str\tstart:int\tend:int\tname:str\tscore:float\tstrand:str"

    for rec in NamedTupleImporter(src,descr=bed6_default,parse_comments=True,default_cast='str',keep_line=True,**kwargs):
        yield rec
    
