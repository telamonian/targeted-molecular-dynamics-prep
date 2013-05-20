import numpy as np
from ..data import Data

class Remarks(Data):
    def __init__(self, pdbparse):
        super(Remarks, self).__init__()
        for remark in pdbparse.rules['REMARK'].data:
            remark = Remark(remark)
            try:
                self.__getattribute__(remark.remarknum).append(remark)
            except AttributeError:
                self.AddDataTuple(((remark.remarknum, [remark]),))
        self.GetBioMatrices()
    
    def GetBioMatrices(self):
        self.biomts = []
        mol_no = 9999
        try:
            for i, remark in enumerate(self.Get('350')):
                if remark.tokens[0:1]==['BIOMT1']:
                    if mol_no > int(remark.tokens[1]):
                        j = self.Get('350')[i-1].tokens.index('CHAINS:')
                        if j>=0:
                            chains = [chain.strip(',') for chain in self.Get('350')[i-1].tokens[j+1:]]
                        else:
                            chains = None
                    mol_no = int(remark.tokens[1])
                    self.biomts.append(BioMT(self.Get('350')[i:i+3], chains))
        except AttributeError:
            pass
            
        
class Remark(object):
    def __init__(self, remark):
        self.remarknum = remark.remarkNum
        self.tokens = remark.empty.split()
        
class BioMT(object):
    def __init__(self, biomt_remarks, chains):
        self.mat = np.identity(4)
        self.chains = chains
        translation = np.identity(4)
        for i, remark in enumerate(biomt_remarks):
            self.mat[i,0:3] = [float(entry) for entry in remark.tokens[2:5]]
            translation[i,3] = -float(remark.tokens[5])
        self.mat = self.mat.dot(translation)
    
    def HasChain(self, chainid):
        try:
            return chainid in self.chains
        except TypeError:
            return True