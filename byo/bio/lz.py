import lz4framed as Z
import os,sys
from time import time
import logging

def chunks(fname, chunksize, alt_src):
    if alt_src:
        f = alt_src
    else:
        f = file(fname,'rb')
    chunk = f.read(chunksize)
    while chunk:
        yield chunk
        chunk = f.read(chunksize)

MB = 1024.**2

from time import time

class ChunkCache(object):
    def __init__(self, max_cached_chunks = 1000):
        self.max_cached_chunks = max_cached_chunks
        
        self.cache = {}
        self.accessed = {}
        self.logger = logging.getLogger('ChunkCache({0})'.format(max_cached_chunks))

    def release(self, fraction=.2):
        """
        free the least accessed fraction (20%) of chunks that are cached to make room for new
        """
        if len(self.cache) < (1-fraction) * self.max_cached_chunks:
            return

        all_keys = [(c,k) for k,c in self.accessed.items()]
        to_drop = sorted(all_keys)[:int(fraction*len(self.accessed))]
        
        self.logger.debug("releasing {0} chunks due to cache limit hit".format(len(to_drop)))
        for c,key in to_drop:
            self.cache.pop(key, None)
            self.accessed.pop(key, None)

        assert len(self.cache) == len(self.accessed)

    def free_lzfile(self, lzfile):
        to_drop = [(lz, chunk_i) for lz, chunk_i in self.cache.keys() if lz == lzfile]

        self.logger.debug("releasing {0} chunks due to LZFile.close()".format(len(to_drop)))
        for c,key in to_drop:
            self.cache.pop(key, None)
            self.accessed.pop(key, None)

        assert len(self.cache) == len(self.accessed)

    def get_chunk(self, lzfile, chunk_i):
        key = (lzfile, chunk_i)
        if not key in self.cache:
            if len(self.cache) > self.max_cached_chunks:
                self.release()

            self.cache[key] = lzfile.get_chunk(chunk_i)
            
        self.accessed[key] = time()
        return self.cache[key]
                
    def __getitem__(self, key):
        lzfile, chunk_i = key
        return self.get_chunk(lzfile, chunk_i)

LZ_chunk_cache = ChunkCache()

class LZFile(object):
    
    chunk_cache = {}
    
    def __init__(self, fname, alt_src = None, level=2, chunk_size=10*MB, compress_on_open=False):
        self.logger = logging.getLogger("LZFile({0})".format(fname) )
        self.basename = fname
        self.chunk_size = chunk_size
        self._pos = 0 # used by seek and __iter__
        ind_file = fname + '.lzot'
        lz_file = fname + '.lzoc'
        if not (os.path.exists(lz_file) and os.path.exists(ind_file)):
            if compress_on_open:
                self.logger.info("compressing '{0}' chunk_size={1} level={2}".format(fname, chunk_size, level))
                self.compress_file(fname, chunksize=chunk_size, alt_src=alt_src, level=level)
            else:
                msg ="The file {0} is not LZ4 compressed and 'compress_on_open' was not set.".format(fname)
                self.logger.error(msg)
                raise IOError(msg)
            

        self.load_index(ind_file)
        self.lz_file = file(lz_file,'rb')

    def close(self):
        self.logger.debug('close() called')
        self.flush()
        
    def flush(self):
        self.logger.debug('flush() called, releasing LZ chunks from global cache')
        LZ_chunk_cache.free_lzfile(self)
        
        
    @staticmethod
    def compress_file( fname, chunksize=10*1024*1024, alt_src = None, level=2):
        import tempfile
        tab_file = tempfile.NamedTemporaryFile(mode="w",dir = os.path.dirname(fname),delete=False)
        tab_file_name = fname + '.lzot'
        comp_file_name = fname + '.lzoc'
        comp_file = file(comp_file_name,'wb')
        comp_base = 0
        cum_size = 0
        t0 = time()
        
        tab_file.write('{0}\n'.format(chunksize))

        for chunk in chunks(fname, chunksize, alt_src = alt_src):
            uncomp_size = len(chunk)

            t1 = time()
            comp = Z.compress(chunk, level=level)
            comp_size = len(comp)
            comp_file.write(comp)
            ratio = 100. * float(comp_size) / uncomp_size 
            t2 = time()
            throughput = cum_size / (t2-t0)

            tab_file.write('{0}\n'.format(comp_base))
            comp_base += comp_size
            cum_size += uncomp_size

            logging.debug("compressed {0}MB ({1:.1f}%) in {2:.1f} sec, {3:.2f} MB/s  sec".format(chunksize/MB, ratio, t2-t1, throughput/MB))

        tab_file.write('{0}\n'.format(comp_base))
        tab_file.write('{0}\n'.format(cum_size))

        # make sure everything is on disk
        os.fsync(tab_file)
        tab_file.close()
                
        # this is atomic on POSIX as we have created tmp in the same directory, 
        # therefore same filesystem
        os.rename(tab_file.name,tab_file_name)

        # make it accessible to everyone
        import stat
        os.chmod(tab_file_name, stat.S_IROTH | stat.S_IRGRP | stat.S_IRUSR)
        os.chmod(comp_file_name, stat.S_IROTH | stat.S_IRGRP | stat.S_IRUSR)

    def load_index(self, idxname):
        self.logger.info("loading index '{0}'".format(idxname))
        lines = file(idxname).readlines()
        self.chunk_size = int(lines[0])
        #self.max_cached_chunks = self.max_cached_MB / (self.chunk_size / 1024**2)
        
        self.L = int(lines[-1])
        self.chunk_starts = [int(b) for b in lines[1:-1]]
        self.logger.info("done.")
        
    def get_chunk(self, i):
        self.lz_file.seek(self.chunk_starts[i])
        comp = self.lz_file.read(self.chunk_starts[i+1] - self.chunk_starts[i])
        return Z.decompress(comp)
    
    def get_chunk_cached(self, i):
        return LZ_chunk_cache.get_chunk(self, i)
            
    def seek(self, pos):
        self._pos = pos

    def __getslice__(self, start, end):
        cs = self.chunk_size
        out = []
        for chunk_i in range(start / cs, (end / cs) + 1):

            chunk_start = chunk_i * cs
            
            c_start = max(start - chunk_start, 0)
            c_end = min(end - chunk_start, cs)
            #print "CHUNK_I", chunk_i, c_start, c_end, cs
            
            out.append(self.get_chunk_cached(chunk_i)[c_start:c_end])
            
        return "".join(out)

    def __iter__(self):
        self.logger.debug("iterating over lines")
        lines = [""]
        cs = self.chunk_size
        
        # simulate seek() functionality
        chunk_i = self._pos / cs
        for i in xrange(chunk_i, len(self.chunk_starts) -1):
            chunk = self.get_chunk_cached(i)
            self.logger.debug("iterating over chunk of {0} bytes".format(len(chunk)))
            
            if i == chunk_i:
                # skip the bytes we are not interested in from the first chunk
                chunk_start = chunk_i * cs
                chunk = chunk[self._pos - chunk_start:]
            
            # and line-iterate over the rest
            new_lines = chunk.splitlines(True)
            if not lines[-1].endswith('\n'):
                # last line from previous chunk was incomplete.
                lines[-1] += new_lines.pop(0)
                
            lines.extend(new_lines)
            for l in lines[:-1]:
                #self.logger.debug("yielding line of {0} bytes".format(len(l)))
                yield l
                
            lines = lines[-1:]

        for l in lines:
            yield l
            
