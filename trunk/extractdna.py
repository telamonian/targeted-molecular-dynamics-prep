'''
Created on Apr 5, 2013

@author: tel
'''
import sys
from pdbparse.pdbparse import PDBParse

import pymol
pymol.pymol_argv = ['pymol', '-q'] #c=command line (no gui), q=quiet (no messages)
pymol.finish_launching()
from pymol import cmd

def Dimerize(pdb):
    fpath = pdb
    pdbparse = PDBParse('pdbparse/pdb_spec_scraped.txt', fpath)
    prot = pdbparse.prot
    remarks = pdbparse.remarks
    prot.MakeBioUnit(remarks.biomts)
    return pdbparse

if __name__=="__main__":
    pdbs = ['2vz4',  '1r8e', '3d6y', '3d6z', '3q1m', '3q2y', '3q3d', '3q5p', '3q5s'] # '3q5r',
    dimers = []
    for pdb in pdbs:
        dimer = Dimerize(pdb).prot
        dimer.PymolLoad(biounit = True)
        dimers.append(dimer)
    
    for dimer in dimers[1:]:
        dimers[0].PymolAlign(dimer, sel1 = 'C. A+C AND i. 3-101', sel2 = 'C. A+C AND i. 2-106')
    for dimer in dimers:
        dimer.PymolSave(filename = dimer.filename + '_aligned.pdb')