'''
Created on Aug 5, 2011

@author: tel
'''
from rule import RuleBase, FWMixin, AutoFWMixin, TokenizeMixin
import unittest, parse
STEP2PATH = '../step2_out.pdb'
PKPATH = '../pK.out'
SUMCRGPATH = '../sum_crg.out'
FORT38PATH = '../fort.38'
HEAD3PATH = '../head3.lst'
RUNPATH = '../run.prm'

class SumcrgIgnoreMixin(object):
    
    def Ignore(self, line):
        if line[0:10] == '----------' or line[0:10] == 'Net_Charge' or line[0:10] == 'Protons   ' or line[0:10] == 'Electrons ':
            return True
        else:
            return False

class Test(unittest.TestCase):
    
    def setUp(self):
        self.lines = ['ATOM   1076  C   CYD A0039_000  13.828  11.422  -7.837   1.700       0.550      BK____M000',
                 'ATOM   1077  O   CYD A0039_000  13.966  10.333  -8.420   1.520      -0.550      BK____M000',
                 'ATOM   1078  CB  CYD A0039_001  13.489  10.672  -5.502   1.700       0.000      01O000M000',
                 'ATOM   1079 1HB  CYD A0039_001  14.471  11.025  -5.257   1.200       0.000      01O000M000',
                 'ATOM   1080 2HB  CYD A0039_001  13.562   9.675  -5.888   1.200       0.000      01O000M000']
        step2 = open(STEP2PATH, 'r')
        pk = open(PKPATH, 'r')
        sumcrg = open(SUMCRGPATH, 'r')
        fort38 = open(FORT38PATH, 'r')
        head3 = open(HEAD3PATH, 'r')
        run = open(RUNPATH, 'r')
        self.step2dict = {'data_tuples': [('natom',(6, 11))],
                          'files': [step2],
                          'rtype': 'Step2',
                          'ppc': (0, 4),
                          'ppm': 'ATOM'}
        self.pkdict = {'files': [pk],
                       'rtype': 'Pk'}
        self.sumcrgdict = {'files': [sumcrg],
                           'rtype': 'Sumcrg'}
        self.fort38dict = {'files': [fort38],
                           'rtype': 'Fort38'}
        self.head3dict = {'files': [head3],
                           'rtype': 'Head3'}
        self.rundict = {'files': [run],
                        'rtype': 'Run',
                        'strip': True,
                        'ptok' : '^\(.*\)$',
                        'mtok': -1,
                        'dtok': 0}

    def tearDown(self):
        for fi in self.step2dict['files']:
            fi.close()
        for fi in self.pkdict['files']:
            fi.close()
        for fi in self.sumcrgdict['files']:
            fi.close()
        for fi in self.fort38dict['files']:
            fi.close()
        for fi in self.head3dict['files']:
            fi.close()
        for fi in self.rundict['files']:
            fi.close()

    def testTokenizeRule(self):
        TokenizeRule = type('TokenizeRule', (RuleBase, TokenizeMixin), {})
        self.runrule = TokenizeRule(**self.rundict)
        self.assertEqual(self.runrule.data[0].EPSILON_PROT, '4.0')

    def testFWRule(self):
        FWRule = type('FWRule', (RuleBase, FWMixin), {})
        self.step2rule = FWRule(**self.step2dict)
        self.assertEqual(int(self.step2rule.data[1000].natom), 1001)
        
    def testAutoFWRule1(self):
        AutoFWRule = type('AutoFWRule', (RuleBase, AutoFWMixin), {})
        self.pkrule = AutoFWRule(**self.pkdict)
        self.assertEqual(self.pkrule.data[10].pH, 'LYS+I0025_')
        self.assertEqual(float(self.pkrule.data[10].dsol), 0.75)
    
    def testAutoFWRule2(self):
        AutoFWRule = type('AutoFWRule', (RuleBase, AutoFWMixin, SumcrgIgnoreMixin), {})
        self.sumcrgrule = AutoFWRule(**self.sumcrgdict)
        self.assertEqual(self.sumcrgrule.data[10].pH, 'LYS+I0025_')
        self.assertEqual(float(self.sumcrgrule.data[10].Get('7')), 1.0)
        
    def testAutoFWRule3(self):
        AutoFWRule = type('AutoFWRule', (RuleBase, AutoFWMixin), {})
        self.fort38rule = AutoFWRule(**self.fort38dict)
        self.assertEqual(self.fort38rule.data[10].ph, 'LYS+1I0005_002')
        self.assertEqual(float(self.fort38rule.data[10].Get('7.0')), 1.0)
        
    def testAutoFWRule4(self):
        AutoFWRule = type('AutoFWRule', (RuleBase, AutoFWMixin), {})
        self.head3rule = AutoFWRule(**self.head3dict)
        self.assertEqual(self.head3rule.data[10].CONFORMER, 'LYS+1I0005_002')
        self.assertEqual(float(self.head3rule.data[10].Get('dsolv')), 0.624)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
