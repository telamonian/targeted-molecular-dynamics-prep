'''
Created on Apr 15, 2012

@author: tel
'''
import re
from data.data import Data
from data.rule import RuleBase, TokenizeMixin, FWMixin, AutoFWMixin, OppMixin, CSVMixin

class PostNameMixin(object):
    '''
    This class implements post processing to munch a standardish mcce id tag into
    resn, res, etc.
    requries:
    self.name_attr  the attribute used by that particular data type to store said MCCE string
    self.nsub_tups  a list of tuples in the form ('name_sub_part', (start index, end index))
    '''
    def PostProcess(self, datum):
        #this mess gives the correct behavior for data from pK.out files
        #since the name column can be called either 'pH' or 'eH'
        if isinstance(self.name_attr, tuple):
            for nattr in self.name_attr:
                try:
                    name = datum.PopDataTuple(nattr)
                    datum.AddDataTuple([('typ', nattr)])
                    break
                except AttributeError:
                    pass
        else:        
            name = datum.PopDataTuple(self.name_attr)
        if hasattr(self, 'nsub_tups'):
            for subname, indexes in self.nsub_tups:
                dtemp = (subname, self.Slice(name, indexes))
                datum.AddDataTuple([dtemp])
        if hasattr(self, 'name_pat'):
            datum.AddDataTuple(re.search(self.name_pat, name).groupdict().items())
        if hasattr(self, 'file_pat'):
            fname = datum.PopDataTuple(self.fname)
            datum.AddDataTuple((self.file_pat[0], re.search(self.file_pat[1], fname).group(0)))

class MCCECmpMixin(object):
    '''
    This class implements __eq__ and __ne__ methods for everything lower
    than runs in the MCCE data object heirarchy
    '''
    def __eq__(self, other):
        foundone = False
        for attr in ('chainid', 'res'):
            try:
                if self.__getattribute__(attr)==other.__getattribute__(attr):
                    foundone = True
                    continue
                else:
                    return False
            except AttributeError:
                continue
        for attr in ('resn', 'confn', 'atomn'):
            try:
                if int(self.__getattribute__(attr))==int(other.__getattribute__(attr)):
                    foundone = True
                    continue
                else:
                    return False
            except AttributeError:
                continue
        return foundone

    def __ne__(self, other):
        return not self.__eq__(other)    

class MCCEData(Data, MCCECmpMixin):
    pass

class MCCERule(RuleBase):
    def __init__(self, **kwargs):
        super(MCCERule, self).__init__(DClass=MCCEData, **kwargs)

class RunRule(MCCERule, TokenizeMixin):
    rtype = 'run.prm'
    def __init__(self, fi):
        self.rundict = {'files': [fi],
                        'rtype': 'run.prm',
                        'strip': True,
                        'ptok' : '^\(.*\)$',
                        'mtok': -1,
                        'dtok': 0}
        super(RunRule, self).__init__(**self.rundict)

class SeqadvRule(MCCERule, FWMixin):
    rtype = 'seqadv'
    def __init__(self, fi):        
        self.seqadvdict = {'data_tuples': [('idcode',(7, 11)), 
                                          ('res',(12, 15)), 
                                          ('chainid',(16, 17)),
                                          ('resn',(18, 22)),
                                          ('icode',(22, 23)),
                                          ('database',(24, 28)),
                                          ('dbidcode',(29, 38)),
                                          ('dbres',(39, 42)),
                                          ('dbseq',(43, 48)),
                                          ('conflict',(49, 70))],
                          'files': [fi],
                          'rtype': 'seqadv',
                          'ppc': (0, 6),
                          'ppm': 'SEQADV'}
        super(SeqadvRule, self).__init__(**self.seqadvdict)

class Step2Rule(MCCERule, FWMixin):
    rtype = 'step2_out.pdb'
    def __init__(self, fi):        
        self.step2dict = {'data_tuples': [('atomn',(6, 11)), 
                                          ('atom',(12, 16)), 
                                          ('altloc',(16, 17)),
                                          ('res',(17, 20)),
                                          ('chainid',(21, 22)),
                                          ('resn',(22, 26)),
                                          ('confn',(27, 30)),
                                          ('x',(30, 38)),
                                          ('y',(38, 46)),
                                          ('z',(46, 54)),
                                          ('rad',(54, 62)),
                                          ('crg',(62, 74)),
                                          ('hist',(74, 90))],
                          'files': [fi],
                          'rtype': 'step2_out.pdb',
                          'ppc': (0, 4),
                          'ppm': 'ATOM'}
        super(Step2Rule, self).__init__(**self.step2dict)

class Head1Rule(MCCERule, FWMixin):
    rtype = 'head1.lst'
    def __init__(self, fi):        
        self.head1dict = {'data_tuples': [('res',(0, 3)), 
                                          ('chainid',(4, 5)), 
                                          ('resn',(5, 9))],
                          'files': [fi],
                          'rtype': 'head1.lst',
                          'ppc': (0, 0),
                          'ppm': ''}
        super(Head1Rule, self).__init__(**self.head1dict)


class Head2Rule(MCCERule, AutoFWMixin, PostNameMixin):
    rtype = 'head2.lst'
    name_attr = 'CONFORMER'
    nsub_tups = [('res',(0,3)),
                 ('conft',(3,5)),
                 ('chainid',(5,6)),
                 ('resn',(6,10)),
                 ('confn',(11,14))]
    def __init__(self, fi):  
        self.head2dict = {'files': [fi],
                          'rtype': 'head2.lst'}
        super(Head2Rule, self).__init__(**self.head2dict)

class Head3Rule(MCCERule, AutoFWMixin, PostNameMixin):
    rtype = 'head3.lst'
    name_attr = 'CONFORMER'
    nsub_tups = [('res',(0,3)),
                 ('conft',(3,5)),
                 ('chainid',(5,6)),
                 ('resn',(6,10)),
                 ('confn',(11,14))]
    def __init__(self, fi):  
        self.head3dict = {'files': [fi],
                          'rtype': 'head3.lst'}
        super(Head3Rule, self).__init__(**self.head3dict)
        
class OppRule(MCCERule, OppMixin):
    rtype = 'energies.opp'
    def __init__(self, fi):  
        self.oppdict = {'files': [fi],
                       'rtype': 'energies.opp'}
        super(OppRule, self).__init__(**self.oppdict)

class FortRule(MCCERule, AutoFWMixin, PostNameMixin):
    rtype = 'fort.38'
    name_attr = ('ph','eh')
    nsub_tups = [('res',(0,3)),
                 ('conft',(3,5)),
                 ('chainid',(5,6)),
                 ('resn',(6,10)),
                 ('confn',(11,14))]
    def __init__(self, fi):  
        self.fortdict = {'files': [fi],
                         'rtype': 'fort.38'}
        super(FortRule, self).__init__(**self.fortdict)

class PKRule(MCCERule, AutoFWMixin, PostNameMixin):
    rtype = 'pK.out'
    name_attr = ('pH','Eh')
    nsub_tups = [('res',(0,3)),
                 ('charge',(3,4)),
                 ('chainid',(4,5)),
                 ('resn',(5,9))]
    def __init__(self, fi):  
        self.pkdict = {'files': [fi],
                       'rtype': 'pK.out'}
        super(PKRule, self).__init__(**self.pkdict)
        
    '''overide the regular self.Process() in order to deal with the 'titration too sharp' lines'''
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
                line = re.sub('titration too sharp', 'titr               ', line)
                line = re.sub('titration curve too sharp', 'titr                     ', line)
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

class CSARule(RuleBase, CSVMixin):
    rtype = 'csa.dat'
    def __init__(self, fi):  
        self.csadict = {'files': [fi],
                       'rtype': 'csa.dat',
                       'headers': ('pdb','site_number','res','chainid','resn','chemical_function','evidence_type','literature_entry')}
        super(CSARule, self).__init__(**self.csadict)

class ConsurfRule(RuleBase, CSVMixin, PostNameMixin):
    rtype = 'consurf.grades'
    name_attr = '3latom'
    name_pat = '(?P<res>[A-Z]{3})(?P<resn>[0-9]+):(?P<chainid>[\w])'
    addfname = True
    file_pat = ('pdb', '/?(\w+)/\w/consurf.grades')
    def __init__(self, fi):  
        self.consurfdict = {'files': [fi],
                       'rtype': 'consurf.grades',
                       'headers': ('pos', 'seq', '3latom', 'score', 'color', 'confidence_interval', 'confidence_interval_colors', 'msa_data', 'residue_variety'),
                       'mtok': 0,
                       'ptok': '^[0-9]+$',
                       'csv_args': {'delimiter': '\t',
                                    'skipinitialspace': True}}
        super(ConsurfRule, self).__init__(**self.consurfdict)