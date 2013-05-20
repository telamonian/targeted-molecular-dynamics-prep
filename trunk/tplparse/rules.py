'''
Created on Mar 26, 2013

@author: tel
'''
import re
from data.data import Data
from data.rule import RuleBase, TokenizeMixin, FWMixin

class BildRule(RuleBase, FWMixin):
    rtype = 'bild'
    def __init__(self, fi):        
        self.bilddict = {'data_tuples': [('atom1',(7, 11)), 
                                          ('atom2',(12, 15)), 
                                          ('atom3',(16, 17)),
                                          ('atom4',(18, 22)),
                                          ('bond12',(22, 23)),
                                          ('angle123',(24, 28)),
                                          ('phi1234',(29, 38)),
                                          ('angle234',(39, 42)),
                                          ('bond34',(43, 48))],
                          'files': [fi],
                          'rtype': 'bild',
                          'ppc': (0, 4),
                          'ppm': 'BILD'}
        super(BildRule, self).__init__(**self.bilddict)
        
if __n