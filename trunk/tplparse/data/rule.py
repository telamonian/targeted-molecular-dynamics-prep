'''
Created on Jul 17, 2011

@author: tel
'''
import re, shlex, csv
import zlib, struct #imports for OppRule
from scipy.sparse import coo_matrix
from data import Data
from collections import namedtuple

class RuleBase(Data):
    '''
    the base rule has no default behavior. simple fixed width behavior. 
    Interface: BaseRule(data_tuples = [(data_name, data_columns)], 
                        files = list of files to data, 
                        lines = any extra lines to include,
                        rtype = type of data produced, defined by self.PreProcess, .Process, [and .Ignore] 
                                (which are in turn defined by mixins) (First Letter Capitalized))
                        DClass = either the Data class, or a subclass, for holding the actual data
                                 (the idea will be to subclass to include sqlalchemy tables)
    '''
    def __init__(self, data_tuples = None, files = None, lines = None, rtype = None, DClass = Data, **kwargs):
        if data_tuples == None:
            data_tuples = []
        if files == None:
            files = []
        if lines == None:
            lines = []
        if rtype == None:
            rtype = 'Base'
        self.DClass = DClass
        self.files = files
        self.lines = lines
        self.rtype = rtype
        super(RuleBase, self).__init__(data_tuples = data_tuples, **kwargs)
        self.PreProcesses()
        if hasattr(self, 'Ignore'):
            self.Ignores()
        self.Processes()
        if hasattr(self, 'PostProcess'):
            self.PostProcesses()

    """deals with adding lines to self.lines. Can handle individual strings or iters of strings"""
    def AddLine(self, line):
        #the condition is True if line has an __iter__ attr, False otherwise
        if getattr(line, '__iter__', False):
            self.lines.append(line)
        else:
            self.lines.append([line])
    
    def PreProcesses(self):
        for lines in self.PreProcess():
            self.AddLine(lines)
    
    '''removes unwanted lines from self.lines
    the reverses and weird indices are due to the fact that it's deleting lines from a list that it's also iterating over'''
    def Ignores(self):
        len1 = len(self.lines) - 1
        for i1, block in enumerate(reversed(self.lines)):
            len2 = len(block) - 1
            for i2, line in enumerate(reversed(block)):
                if self.Ignore(line):
                    block.pop(len2 - i2)
            if block == []:
                self.lines.pop(len1 -i1)
    
    def Processes(self):
        for dtups in self.Process():
            datum = self.DClass(dtype=self.rtype, data_tuples=dtups)
            self.AddData(datum)
    
    def PostProcesses(self):
        for datum in self.data:
            self.PostProcess(datum)
    
    @staticmethod
    def Slice(line, col_tup, strip = True):
        if strip:
            return line[col_tup[0]:col_tup[1]].strip()
        else:
            return line[col_tup[0]:col_tup[1]]

class FWMixin(object):
    '''
    adds simple fixed width behavior to RuleBase. 
    Interface: BaseRule(data_tuples = [(data_name, data_columns)], 
                        files = list of files to data, 
                        lines = any extra lines to include,
                        rtype = name of data to be produced (First Letter Capitalized),
                        ppc = pre-process columns. tuple with slices indices to check for ppm
                        ppm = pre-process match. if found within ppc, the line is included (i.e. ATOM for pdb coords))
    '''
    '''fixed width parser behavior. self.ppc = pre-process columns (a tuple), self.ppm = pre-process match'''
    def PreProcess(self):
        for fi in self.files:
            for line in fi:
                if self.Slice(line, self.ppc) == self.ppm:
                    yield line
    
    '''fixed width process behavior.'''
    def Process(self):
        dtups = []
        #line[0] is needed to 'liberate' the strings that actually represent the file lines from the list of lists that is self.lines
        for line in (line[0] for line in self.lines):
            dtup = []
            for key, value in self.GenDataTuples():
                dtemp = (key, self.Slice(line, value))
                dtup.append(dtemp)
            dtups.append(dtup)
        return dtups
    
class AutoFWMixin(object):
    '''
    Adds automatic fixed width behavior to RuleBase. Can only deal with one file at a time. 
    Interface: BaseRule(data_tuples = [(data_name, data_columns)], 
                        files = list of files to data, 
                        lines = any extra lines to include,
                        rtype = name of data to be produced (First Letter Capitalized),
                        ppc (optional) = pre-process columns. tuple with slices indices to check for ppm
                        ppm (optional) = pre-process match. if found within ppc, the line is included (i.e. ATOM for pdb coords))
                        ignore (optional) = list of functions that take a line and return true if it is to be ignored
    '''
    '''fixed width parser behavior. self.ppc = pre-process columns (a tuple), self.ppm = pre-process match'''
    def PreProcess(self):
        if not hasattr(self, 'Ignore'):
            self.Ignore = lambda x: False
        header = self.files[0].readline()
        self.files[0].seek(0)
        always_space = [1]*len(header)
        always_transition = [1]*(len(header))
        for line in self.files[0]:
            if self.Ignore(line):
                continue
            else:
                lastchar = ' '
                for i, (space_spot, trans_spot, char) in enumerate(zip(always_space, always_transition, line)):
                    always_space[i] = char.isspace() & space_spot
                    always_transition[i] = (char.isspace() | lastchar.isspace()) & trans_spot
                    lastchar = char
        data_tuples = []
        start_slice = None
        for i, (space_spot, trans_spot) in enumerate(zip(always_space, always_transition)):
            if start_slice != None:
                if trans_spot==1:
                    btup = (start_slice, i)
                    data_tuples.append((self.Slice(header, btup), btup))
                    if space_spot!=1:
                        start_slice = i
                    else:
                        start_slice = None
                else:
                    continue
            elif trans_spot==1 and space_spot!=1:
                start_slice=i
            else:
                continue
        self.AddDataTuple(data_tuples)
        self.files[0].seek(0)
        self.files[0].readline()
        for line in self.files[0]:
            yield line
        
    '''fixed width process behavior.'''
    def Process(self):
        dtups = []
        for line in (line[0] for line in self.lines):
            dtup = []
            for key, value in self.GenDataTuples():
                dtemp = (key, self.Slice(line, value))
                dtup.append(dtemp)
            dtups.append(dtup)
        return dtups

class CSVMixin(object):
    '''
    Adds comma separated value behavior to RuleBase. Can only deal with one file at a time. 
    Interface: BaseRule(data_tuples = [(data_name, data_columns)], 
                        files = list of files to data, 
                        lines = any extra lines to include,
                        rtype = name of data to be produced (First Letter Capitalized),
                        csv_args (optional) = keyword arguments to be passed to the csv.reader object
                        headers (optional) = list of headers. Used for files that lack header lines. To overide existing file header lines,
                                             specify self.headers and an appropriate self.Ignore function.
                        nohead (optional) = If this is True, then the first line of the file is assumed not to be a
                                               header, and is included where otherwise it would be ignored
                        Ignore (optional) = list of functions that take a line and return true if it is to be ignored
                        ptok (optional) = token pattern. If the regex pattern in ptok is found in the match token (mtok), the line is claimed
                        mtok (optional) = match token. Index of token to check for ptok.
                        addfname (optional) = adds the name of the originating file to every datum as an fname property
    '''
    def PreProcess(self):
        if not hasattr(self, 'csv_args'):
            self.csv_args = {'delimiter':','}
        if not hasattr(self, 'headers'):
            self.headers = csv.reader([self.files[0].readline()], **self.csv_args).next()
        self.files[0].seek(0)
        if hasattr(self, 'nohead'):
            if self.nohead!=True:
                self.files[0].readline()
        else:
            self.files[0].readline()
        if hasattr(self, 'addfname'):
            if self.addfname==True:
                for line, fname in ((line, fi.name) for fi in self.files for line in fi):
                    yield '%s%s%s' % (line,self.csv_args['delimiter'],fname)
            else:
                for line in (line for fi in self.files for line in fi):
                    yield line
        else:
            for line in (line for fi in self.files for line in fi):
                yield line

    def Process(self):
        dtups = []
        for values in csv.reader((line[0] for line in self.lines), **self.csv_args):
            values = filter(None, values)
            if hasattr(self, 'ptok') and hasattr(self, 'mtok'):
                if re.search(self.ptok, values[self.mtok]):
                    dtup = zip(self.headers, values)
                    dtups.append(dtup)
                else:
                    continue
            else:
                dtup = zip(self.headers, values)
                dtups.append(dtup)
        return dtups

class TokenizeMixin(object):
    '''
    adds token based parsing to rule base.
    Interface: BaseRule(data_tuples = [(data_name, data_columns)], 
                        files = list of files to data, 
                        lines = any extra lines to include,
                        rtype = name of data to be produced (First Letter Capitalized),
                        ptok = token pattern. If the regex pattern in ptok is found in the match token (mtok), the line is claimed
                        mtok = match token. Index of token to check for ptok. Becomes the data header
                        dtok = data token. Index of token to take as data in matched line. -1, of course, matches the last token
                        strip = optional logical value. If True, strip out parentheses from start and end of tokens
    '''
    def PreProcess(self):
        for fi in self.files:
            for li in fi:
                try:
                    tokens = shlex.split(li)
                except ValueError:
                    tokens = li.split()
                if tokens:
                    if re.search(self.ptok, tokens[self.mtok]):
                        yield tokens
    
    def Process(self):
        dtups = []
        for tokens in self.lines:
            mtok = tokens[self.mtok]
            dtok = tokens[self.dtok]
            if getattr(self, 'strip', False):
                mtok = mtok.strip('()')
                dtok = dtok.strip('()')
            dtup = (mtok, dtok)
            dtups.append(dtup)
        return [dtups]
    
class OppMixin(object):
    '''
    adds the ability to parse zipped MCCE energies.opp files to rule base.
    probably only useful for the one type of file.
    Interface: BaseRule(data_tuples = [(data_name, data_columns)], 
                        files = list of files to parse (should only contain one energiess.opp file), 
                        lines = any extra lines to include,
                        rtype = name of data to be produced (First Letter Capitalized),
    '''
    CHUNKSIZE = 2**14   #2**24 is ~32 million, so 32 MB
    CONF_HEAD = struct.Struct('15sc4f2i8f12s')
    CONF_HEAD_fields = ('uniqID', 'on', 'occ', 'netcrg',
                        'Em', 'pKa', 'e', 'H',
                        'E_vdw0', 'E_vdw1', 'E_tors', 'E_epol',
                        'E_rxn0', 'E_rxn', 'E_dsolv', 'E_extra',
                        'history')
    PAIRWISE = struct.Struct('4f4s')
    PAIRWISE_fields = ('ele', 'vdw', 'crt', 'ori', 'mark')
    
    @staticmethod
    def fread(size, count, stream):
        '''
        this is gonna be a little different from the traditional C version
        will return a generator which will yield each sucessive block rather than storing them directly into an array
        '''
        for i in xrange(count):
            yield stream.Read(size)
        
    
    class DFile(object):
        '''
        takes a zlib compressed file, reads in self.chunksize compressed bytes at a time, and then hands back
        the equivalent uncompressed string. Also keeps track of the equivalent uncompressed file pointer.
        inherits from zlib.decompressobj so it shares the same methods.
        '''
        def __init__(self, fi, chunksize=2**8, offset=0):
            '''
            offset can be set to the string 'dynamic'
            this assumes that the header goes from the start of the file until the first newline, i.e. an energies.opp style header
            '''
            self.file = fi
            self.chunksize = chunksize
            self.offset = offset
            self.d = zlib.decompressobj()
            self.FirstChunkRead()
        
        def FirstChunkRead(self):
            self.ChunkRead()
            if self.offset:
                if self.offset == 'dynamic':
                    wspace = self.dbuff.find(' ')
                    nline = self.dbuff.find('\n')    #every energies.opp starts with a header with format 'n_conf MCCE_verno'. This helps locate it
                    self.n_conf, self.verno, self.dbuff = (int(self.dbuff[:wspace]), self.dbuff[wspace+1:nline], self.dbuff[nline+1:])
                else:
                    self.header, self.dbuff = (self.dbuff[:self.offset], self.dbuff[self.offset:])
                     
        def ChunkRead(self, tail=''):
            self.dbuff = tail + self.d.decompress(self.file.read(self.chunksize))
        
        def Read(self, byte):
            '''
            consume byte # of bytes from the uncompressed file
            if the end of dbuff is reached before byte # of bytes are read, preserve dbuff tail, read in next chunk, and append tail to front, then keep going
            '''
            while len(self.dbuff) < byte:
                self.ChunkRead(self.dbuff)
            out, self.dbuff = (self.dbuff[:byte], self.dbuff[byte:])
            return out

    class Conf_Head(object):
        def __init__(self, f):
            self.dfile = OppMixin.DFile(f, offset='dynamic')
        
        def Run(self):
            conf_heads = []
            for i, ch in enumerate(OppMixin.fread(OppMixin.CONF_HEAD.size, self.dfile.n_conf, self.dfile)):
                conf_head = OppMixin.PAIRWISE.unpack_from(ch)
                conf_heads.append(conf_head)
            return conf_heads
    
    class Pairwise(object):
        def __init__(self, f):
            #if f is a file object, prep it and consume up to the end of the CONF_HEAD section
            #otherwise, assume that f is a DFile that has already been consumed by a Conf_Head object up to the end of the CONF_HEAD section
            if isinstance(f, file):
                self.dfile = OppMixin.DFile(f, offset='dynamic')
                #fast forward
                self.dfile.Read(self.dfile.n_conf*OppMixin.CONF_HEAD.size)
            else:
                self.dfile = f
        
        def Run(self):
            ei,ej,ev,vi,vj,vv = [],[],[],[],[],[]
            for i in xrange(self.dfile.n_conf):
                for j, pw in enumerate(OppMixin.fread(OppMixin.PAIRWISE.size, self.dfile.n_conf, self.dfile)):
                    epw, vpw = OppMixin.PAIRWISE.unpack_from(pw)[0:2]
                    if epw != 0:
                        ei.append(i+1)
                        ej.append(j+1)
                        ev.append(epw)
                    if vpw != 0:
                        vi.append(i+1)
                        vj.append(j+1)
                        vv.append(vpw)
            epw = coo_matrix((ev,(ei,ej)),shape=(self.dfile.n_conf+1,self.dfile.n_conf+1))
            vpw = coo_matrix((vv,(vi,vj)),shape=(self.dfile.n_conf+1,self.dfile.n_conf+1))
            epw, vpw = epw.tocsr(), vpw.tocsr()
            dtups = [(('epw',epw), 
                      ('vpw',vpw))]
            return dtups
            
    def PreProcess(self):
        return []   #this skips preprocessing so that the complete lines of this file are never stored in memory
                    #for 15,000 conformers reading the whole thing into memory twice takes 8 gb, so we want to avoid that
    
    def Process(self):
        dtups = []
        for fi in self.files:
            pairwise = OppMixin.Pairwise(fi)
            dtups += pairwise.Run()
        return dtups
                
                
class DummyMixin(object):
    '''
    allows passing of data_tuples directly into the parsing object as a list of lists
    Interface: BaseRule(data_tuples = [(data_name, data_columns)], 
                        data_tuple_list = list of lists that become data objects in data attribute
    '''
    def PreProcess(self):
        return []
    
    def Process(self):
        return self.data_tuple_list
  
#    @staticmethod
#    def BreakAt(line, breaks):
#        broken = []
#        broken.append(line[:breaks[0]])
#        for br1, br2 in zip(breaks[:-1], breaks[1:]):
#            broken.append(line[br1:br2])
#        broken.append(line[breaks[-1]:])
#        return broken   
