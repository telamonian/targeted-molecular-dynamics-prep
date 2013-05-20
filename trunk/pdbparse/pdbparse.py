'''
Created on Sep 20, 2012

@author: tel
'''
import re, os
from StringIO import StringIO
from rule import RuleBase, AutoFWMixin, FWMixin
from urllib import urlretrieve

from pdbobjs.prot import Prot
from pdbobjs.remarks import Remarks

LINEWIDTH = 128

class PDBAutoFWMixin(AutoFWMixin):
    '''
    fixed width process behavior.
    modified for PDBParse objects so that the 'unclaimed' space between columns is assigned to the adjacent column to the left.
    this was done (in conjunction with automatically parsing only the header) in order to deal with the white space in pdb spec entries.
    '''
    def Process(self):
        self.Widen()
        return AutoFWMixin.Process(self)
    
    def Widen(self):
        sorted_data_tuples = sorted(self.GetDataTuples(), key=lambda tup: tup[1])
        for i, (key, value) in enumerate(sorted_data_tuples[:-1]):
            value = (value[0], sorted_data_tuples[i+1][1][0])
            self.__setattr__(key, value)
        lastone = sorted_data_tuples[-1]
        self.__setattr__(lastone[0], (lastone[1][0], None))
    
class PDBParse(object):
    class HeaderRule(RuleBase, PDBAutoFWMixin):
        rtype = 'header'
        def __init__(self, fi):  
            self.specdict = {'files': [fi],
                           'rtype': 'header'}
            super(self.__class__, self).__init__(**self.specdict)  
    
    class BlockRule(RuleBase, FWMixin):
        rtype = 'block'
        def __init__(self, blockdict):        
            self.blockdict = blockdict
            super(self.__class__, self).__init__(**self.blockdict)
            
    class PDBRule(RuleBase, FWMixin):
        rtype = 'PDB'
        def __init__(self, pdbpath, pdbdict):
            '''
            assume that pdbs which cannot be opened can be found online
            '''
            try:
                fi = open(pdbpath)
            except IOError:
                try:
                    pdbpath_corrected = '%s.pdb' % '.'.join(pdbpath.split('.')[0:1] + pdbpath.split('.')[1:-1])
                    fi = open(pdbpath_corrected)
                except IOError:
                    urlretrieve(os.path.join('http://www.rcsb.org/pdb/files', pdbpath_corrected), pdbpath_corrected)
                    fi = open(pdbpath_corrected)
            self.pdbdict = pdbdict
            self.pdbdict['files'] = [fi]
            super(self.__class__, self).__init__(**self.pdbdict)
            fi.close()
    
    def __init__(self, specpath, pdbpath):
        self.rules = self.GetRules(specpath, pdbpath)
        self.pdbpath = '%s.pdb' % '.'.join(pdbpath.split('.')[0:1] + pdbpath.split('.')[1:-1])
        self.prot = Prot(self)
        self.remarks = Remarks(self)
        
    def GetRules(self, specpath, pdbpath):
        rules = {}
        self.specf = open(specpath)
        for header, block in zip(*self.PreParse()):
            blockdict = {'files': [block], 'rtype': 'block', 'ppc': (0,0), 'ppm': ''}
            blockdict['data_tuples'] = self.__class__.HeaderRule(header).GetDataTuples()
            blockrule = self.__class__.BlockRule(blockdict)
            pdbdict = {'data_tuples':[]}
            ppc = re.search('^\s*(\d+)\s*-\s*(\d+)', blockrule.data[0].COLUMNS)
            pdbdict['ppc'] = (int(ppc.group(1)) - 1, int(ppc.group(2)))
            ppm = blockrule.data[0].FIELD.strip('"').strip()
            pdbdict['ppm'] = ppm
            pdbdict['rtype'] = ppm
            for datum in blockrule.data[1:]:
                if not datum.COLUMNS.isspace():
                    cols = re.search('^\s*(\d+)\s*(?:-\s*(\d+))?', datum.COLUMNS)
                    if cols:
                        if cols.group(1) and cols.group(2):
                            pdbdict['data_tuples'].append((datum.FIELD.strip('"'), (int(cols.group(1)) - 1, int(cols.group(2)))))
                        else: #elif cols.group(1): don't need the elif here, but it makes the code more readable
                            pdbdict['data_tuples'].append((datum.FIELD.strip('"'), (int(cols.group(1)) - 1, int(cols.group(1)))))
            rules[pdbdict['rtype']] = self.__class__.PDBRule(pdbpath, pdbdict)
        return rules
        
    def PreParse(self):
        headers = []
        blocks = []
        block = ''
        inblock = False
        for line in self.specf:
            if re.match('^-{20}', line):
                continue
            elif inblock==True and (re.match('^\s*\d+\s+', line) or re.match('\s{10}', line)):
                block += line
            elif inblock==True:
                blocks.append(StringIO(block))
                block = ''
                inblock = False
            elif re.match('^COLUMNS', line):
                line = re.sub('DATA TYPE', 'DATA_TYPE', line)
                line = re.sub('DATA  TYPE', 'DATA__TYPE', line)
                headers.append(StringIO(line))
                block = ''
                inblock = True
            else:
                pass
        return headers, blocks