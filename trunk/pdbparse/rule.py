'''
Created on Jul 17, 2011

@author: tel
'''
from data import Data

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
        if hasattr(self, 'Ignore'):
            self.Ignores()
        self.PreProcesses()
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
                ignore = self.Ignore(line)
                if ignore==None:
                    block.pop(len2 - i2)
                else:
                    block[len2 - i2] = ignore
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
                if self.Slice(line, self.ppc) == self.ppm.strip():
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
        header = self.files[0].readline()
        self.files[0].seek(0)
        always_space = [1]*len(header)
        always_transition = [1]*(len(header))
        for line in self.files[0]:
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
