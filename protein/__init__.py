from codons import codon_table

def translate(seq,frame=0):
    started = False
    for i in xrange(frame,len(seq)-frame-4,3):
        cod = seq[i:i+3].upper()
        a = codon_table.get(cod,"X")

        if not started and a == 'M':
            # scanning for start codon
            started = True

        if started:
            # within the ORF
            yield a

        if started and a == '*':
            # translating and encountered stop codon
            break

def translate_aa(seq,frame=0):
    return "".join(translate(seq,frame))

def all_orfs(seq,need_stop=False):
    for frame in range(3):
        started = False
        start = 0
        orf = []

        for i in xrange(frame,len(seq)-frame-3,3):
            cod = seq[i:i+3].upper()
            a = codon_table.get(cod,"X")

            if (not started) and a == 'M':
                # scanning for start codon
                started = True
                start = i
                orf = []
                #print "* starting at %d" % i,orf

            if started:
                # within the ORF
                orf.append(a)
                #print "* elongating",orf

            if started and a == '*':
                # translating and encountered stop codon
                #print "* terminating",orf
                yield (start,i+3,"".join(orf))
                started = False
                orf = []

        if orf and not need_stop:
            yield (start,i+3,"".join(orf))

def find_ORF(seq,len_thresh=30):
    """
    find all open-reading frames, sort by length (if above threshold) and return the longest one.
    """
    orfs = sorted([ (len(orf),orf,start,end) for (start,end,orf) in all_orfs(seq) if len(orf) >= len_thresh],reverse=True)
    if not orfs:
        return "",-1,-1
    else:
        L,ORF,start,end = orfs[0]
        return ORF,start,end

def trypsinize(aa,minlen=6,undigested=2):
    """
    in-silico digest the amino acid sequence with trypsin, cutting after R or K.
    can allow for arbitrary number of undigested residues to simulate partial
    digestion (default=up to 2).
    """
    import re
    bps = sorted(set([-1,] + [M.span()[0] for M in re.finditer(r'(R|K|r|k)',aa)] + [len(aa)-1,]))

    from itertools import tee, izip
    def n_wise(iterable,n):
        "s -> (s0,sn), (s1,s1+n), (s2, s2+n), ..."
        a, b = tee(iterable)
        for i in range(n):
            next(b, None)
        return izip(a, b)
    
    for j in range(1,undigested+2):
        #print "undigested",j
        for bp1,bp2 in n_wise(bps,j):
            pep = aa[bp1+1:bp2+1]
            if len(pep.replace('*','')) >= minlen:
                yield pep
