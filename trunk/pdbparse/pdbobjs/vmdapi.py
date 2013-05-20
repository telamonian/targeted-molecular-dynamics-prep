'''
Created on Apr 29, 2013

@author: tel
'''
import os
from subprocess import call

SOLVPAD = 5.0
SALTCON = 0.15
CHARMM = 27

PROTRES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
NARES = ['A','C','G','U','DA','DC','DG','DT']
PURINES = ['A','G','DA','DG']
PYRIMIDINES = ['C','U','DC','DT']
WATRES = ['HOH']

VMDPATH = '/Applications/VMD 1.9.1.app/Contents/vmd/vmd_MACOSXX86'

class VMDProt(object):
    def Savepqn(self, dirpath, biounit=False, target=False):
        if CHARMM==27:
            topfiles = ['topology top_all27_prot_na.rtf']
        elif CHARMM==36:
            topfiles = ['topology top_all36_prot.rtf', 
                    'topology top_all36_na.rtf', 
                    'topology toppar_water_ions.str']
        pqnboilerplate = ['package require psfgen',
                    'package require autoionize',
                    'package require solvate']+topfiles+['', 
                    'pdbalias residue DG GUA ', 
                    'pdbalias residue DC CYT ', 
                    'pdbalias residue DA ADE ', 
                    'pdbalias residue DT THY', 
                    'pdbalias residue HIS HSE', 
                    'pdbalias atom ILE CD1 CD', 
                    'pdbalias residue HOH TIP3',
                    'pdbalias atom HOH O OH2', '']
        pqnendplate = ['guesscoord',
                    'writepdb %s_vmd.pdb' % self.name,
                    'writepsf %s_vmd.psf' % self.name,
                    'solvate %s_vmd.psf %s_vmd.pdb -t %.3f -o %s_vmd_wat' % (self.name, self.name, SOLVPAD, self.name),
                    'autoionize -psf %s_vmd_wat.psf -pdb %s_vmd_wat.pdb -sc %.3f -cation SOD -anion CLA -o %s_vmd_wat_ion' % (self.name, self.name, SALTCON, self.name),
                    'quit']
        patches = []
        with open(os.path.join(dirpath,self.name+'.pqn'), 'w') as f:
            for line in pqnboilerplate:
                f.write('%s\n' % line)
            for chain in self.IterChains(biounit, target):
                seglines = ['segment %s {' % chain.chainid, 'pdb %s_%s.pdb' % (self.name.split('_')[0],chain.chainid)]
                resiter = chain.IterRess()
                try:
                    res0 = resiter.next()
                    if res0.resname in NARES:
                        seglines.append('first 5TER')
                        seglines.append('last 3TER')
                        if CHARMM==27:
                            if res0.resname in PYRIMIDINES:
                                patches.append('patch DEO1 %s:%s' % (chain.chainid, res0.resseq))
                            else:
                                patches.append('patch DEO2 %s:%s' % (chain.chainid, res0.resseq))
                            for res in resiter:
                                if res.resname in PYRIMIDINES:
                                    patches.append('patch DEO1 %s:%s' % (chain.chainid, res.resseq))
                                else:
                                    patches.append('patch DEO2 %s:%s' % (chain.chainid, res.resseq))
                        elif CHARMM==36:
                            patches.append('patch DEO5 %s:%s' % (chain.chainid, res0.resseq))
                            for res in resiter:
                                patches.append('patch DEOX %s:%s' % (chain.chainid, res.resseq))
                    elif res0.resname in WATRES:
                        seglines.append('auto none')
                except StopIteration:
                    pass
                seglines += ['}', '']   #the empty string creates a linebreak between segment commands
                seglines.append('coordpdb %s_%s.pdb %s' % (self.name.split('_')[0], chain.chainid, chain.chainid))
                for line in seglines:
                    f.write('%s\n' % line)
            for patch in patches:
                f.write('%s\n' % patch)
            for line in pqnendplate:
                f.write('%s\n' % line)
    
    def Psfgen(self, dirpath, biounit=False, target=False, overwrite=True):
        fpath = os.path.join(dirpath,'%s.pqn' % self.name)
        if not os.path.isfile(fpath) or overwrite:
            self.Savepqn(dirpath, biounit, target)
        os.chdir(dirpath)
        cmdtok = [VMDPATH,'-dispdev', 'text', '-e', '%s.pqn' % self.name]
        call(cmdtok)
        os.chdir('..')