'''
Created on Mar 18, 2013

@author: tel
'''
import os, time, sys
from pdbparse.pdbparse import PDBParse

PDBSPEC = 'pdbparse/pdb_spec_scraped.txt'

def InitPymol():
    '''
    setup pymol so that it works right and plays nice with this program
    '''
    import pymol
    pymol.pymol_argv = ['pymol', '-q'] #c=command line (no gui), q=quiet (no messages)
    pymol.finish_launching()
    from pymol import cmd
    for switchon in ['retain_order','pdb_no_end_record','cartoon_cylindrical_helices']: #'pdb_retain_ids',
        cmd.set(switchon, 1)
    for switchoff in ['pdb_use_ter_records']:
        cmd.set(switchoff, 0)
    cmd.show_as('cartoon')

def LoadBiounit(pdb):
    '''
    pdb = name of pdb file, minus the .pdb filename extension, to load into a Prot object and make the biounit of
    if file is in working directory, use that
    otherwise, download the pdb from pdb.org to working directory and use that
    EXAMPLE: LoadBiounit('3d6z')
    '''
    fpath = pdb
    pdbparse = PDBParse(PDBSPEC, fpath)
    prot = pdbparse.prot
    remarks = pdbparse.remarks
    prot.MakeBioUnit(remarks.biomts)
    prot.SepHet(biounit=True)
    prot.RemoveNamelist(['GOL','RHQ','KAN','IMD','P4P','ACH','TAC','PUY','ET','M4A','BER'], biounit=True)
    prot.RemoveFarWater(cutoff=4, biounit=True)
    prot.Center(biounit=True)
    prot.Orient(biounit=True)
    prot.RenumberNegs(biounit=True)
    return prot

def MakeAlign(p1, p2, sel1, sel2, biounit=True):
    '''
    p1 = Prot object used as the mobile structure in the alignment, 
    p2 = Prot object used as the stationary structure in the alignment, 
    sel1 = pymol select string for mobile, 
    sel2 = pymol select string for stationary
    EXAMPLE: MakeAlign(p1=LoadBiounit('2vz4'), p2=LoadBiounit('3d6z'), sel1='C. A+C AND i. 2-106', sel2='C. A+C AND i. 3-101')
    '''
    for prot in (p1, p2):
        prot.PymolLoad(biounit)
    p2.PymolAlign(p1, sel1, sel2)
    for prot in (p1, p2):
        prot.PymolSave(filename=os.path.join(prot.name, prot.name + '_aligned.pdb'))

def MakePSF(pdb):
    alignedprot = PDBParse(PDBSPEC, os.path.join(pdb, pdb + '_aligned.pdb')).prot
    alignedprot.SepHet() 
    #alignmer.Transform(alignmer.GetOrientTransmat())
    alignedprot.SaveChains(os.path.join(alignedprot.name.split('_')[0],alignedprot.name.split('_')[0]))
    alignedprot.Psfgen(dirpath=pdb)

def MakeTargets(startpdb, endpdb):
    vmdstartprot = PDBParse(PDBSPEC, os.path.join(startpdb, startpdb + '_aligned_vmd.pdb')).prot
    vmdstartprot.SepHet()
    alignedendprot = PDBParse(PDBSPEC, os.path.join(endpdb, endpdb + '_aligned.pdb')).prot
    vmdstartprot.MakeTargetPDB(alignedendprot, outpath=os.path.join(startpdb, startpdb + '_aligned_vmd_target.pdb'))
    alignedendprot.MergeSolvTarget(os.path.join(startpdb, startpdb + '_aligned_vmd_target.pdb'), 
                                   os.path.join(startpdb, startpdb + '_aligned_vmd_wat_ion.pdb'))
    #for illustration purposes
    vmdstartprot.MakeTargetPDB(alignedendprot, os.path.join(startpdb, startpdb + '_aligned_vmd_target_ca.pdb'), ca=True)
    vmdstartprot.Save(os.path.join(startpdb, startpdb + '_aligned_vmd_ca.pdb'), ca=True)

def GetSolvDims(pdb):
    '''
    Get total system dimensions (vec of min, vec of max, difference of the two, average of the two)
    of the solvated, salted system
    '''
    prot = PDBParse(PDBSPEC, os.path.join(pdb, pdb + '_aligned_vmd_wat_ion.pdb')).prot
    return prot.GetDims()
    
def PrepForTMD(startpdb, endpdb, startsel, endsel):
    '''
    the main loop for the program
    prepares the starting structure for NAMD in the expected way, producing a solvated, salted pdb and psf of the biounit
    prepares the end structure as a NAMD-appropriate target pdb
    should work right whether the end structure is identical to the start, or if the end structure is a fragment of a homolog of the start, or anywhere inbetween
    emphasis on the should, so double check the NAMD input files before starting a proper simulation
    startpdb = starting structure pdb filename for TMD. Leave '.pdb' off. If it is not in the working directory, it will be downloaded from pdb.org
    endpdb = end structure pdb filename for TMD. Leave '.pdb' off. If it is not in the working directory, it will be downloaded from pdb.org
    startsel = pymol select string for start structure. Used during the alignment of the end structure to the start structure
    endsel = pymol select string for end structure. Used during the alignment of the end structure to the start structure
    '''
    for pdb in (startpdb, endpdb):
        try:
            os.makedirs(pdb)
        except OSError:
            pass
    now = time.time()
    MakeAlign(p1=LoadBiounit(endpdb), p2=LoadBiounit(startpdb), sel1=endsel, sel2=startsel)
    print 'Alignment of target to start took %.3f seconds' % (time.time() - now)
    now = time.time()
    MakePSF(startpdb)
    print 'Creation of salted, solvated .pdb/.psf of start structure took %.3f seconds' % (time.time() - now)
    now = time.time()
    MakeTargets(startpdb, endpdb)
    print 'Creation of target .pdb took %.3f seconds' % (time.time() - now)
    now = time.time()
    print GetSolvDims(startpdb)
    print 'Calculation of solvated, salted dimensions took %.3f seconds' % (time.time() - now)
    
if __name__=="__main__":
    startpdb, endpdb, startsel, endsel = sys.argv[1:5]
    InitPymol()
    PrepForTMD(startpdb, endpdb, startsel, endsel)
#    pdbs = ['2vz4','3d6z']# '3d6y', '3q1m', '3q2y', '3q3d', '3q5p', '3q5s','1r8e','3q5r'] #['3cyt','1csv']
#    dimers = []
#    for pdb in pdbs:
#        try:
#            os.makedirs(pdb)
#        except OSError:
#            pass
#        dimer = Dimerize(pdb).prot
#        dimer.PymolLoad(biounit = True)
#        dimers.append(dimer)
#    now = time.time()
#    for dimer in dimers[1:]:
#        dimer.PymolAlign(dimers[0], sel1 = 'C. A+C AND i. 2-106', sel2 = 'C. A+C AND i. 3-101')
#    for dimer in dimers:
#        dimer.PymolSave(filename = os.path.join(dimer.filename,dimer.filename + '_aligned.pdb'))
#    alignmers = []
#    pdb_aligned = [os.path.join(pdb,'%s_aligned.pdb' % pdb) for pdb in pdbs]
#    for pdb in pdb_aligned:
#        alignmers.append(PDBParse(PDBSPEC, pdb).prot)
#    for alignmer in alignmers:
#        alignmer.SepHet()
#        #alignmer.Transform(alignmer.GetOrientTransmat())
#        alignmer.SaveChains(os.path.join(alignmer.filename.split('_')[0],alignmer.filename.split('_')[0]))
#    alignmers[1].MakeTargetPDB(alignmers[0], '3d6z/3d6z_target.pdb')
#    alignmers[1].Psfgen(dirpath='3d6z')
#    alignmers[1].MergeSolvTarget('3d6z/3d6z_target.pdb', '3d6z/3d6z_aligned_vmd_wat_ion.pdb')
#    #for illustration purposes
#    alignmers[1].MakeTargetPDB(alignmers[0], '3d6z/3d6z_target_ca.pdb', ca=True)
#    alignmers[1].Save('3d6z/3d6z_ca.pdb', ca=True)
#    print time.time() - now

#    for dimer in dimers[1:]:
#        result = dimers[0].PairMin(dimer, srange_tups=(('A',(2,102)), ('C',(2,102))), orange_tups=(('A',(2,107)), ('C',(2,107))), biounit=True, ca=True)
#        #result = dimers[0].PairMin(dimer, srange_tups=(('A',(2,102)), ('C',(2,102)), ('B',(-10,10)), ('D',(-10,10))), orange_tups=(('A',(6,107)), ('C',(6,107)), ('B',(-10,10)), ('D',(-10,10))), biounit=True, ca=True)
#        #result = dimers[0].PairMin(dimer, srange_tups=(('O',(1,100)),), orange_tups=(('A',(1,100)),), biounit=True, ca=True)
#        print result
#        angvecpnt = result[1]
#    dimers[0].RotateSplit(angvecpnt[0], angvecpnt[1:4], angvecpnt[4:7], biounit=True, ca=True)
#    dimers[0].Rotate(2.775980963459502, [0.12353733144956269, 0.3747743641710698, 0.8042689819549111], [12.009845427786509, -18.127282690241174, 6.8810274778421], biounit=True)
#    for dimer in dimers:
#        dimer.Save(fpath=dimer.filename + '_maxlign.pdb', biounit=True)

    #dimers[0].Objective(dimers[1], [0,0,0,0,0,0], biounit=False, ca=True)
    #transrot = [13.93815383, -23.30498079, 9.44130034] + [2.5,5,2.5]
#    dimers[0].RotateBiounit(*dimers[0].MinVecSplit(dimers[1], biounit=True))
    #dimers[0].RotateBiounit(transrot[0:3], transrot[3:6])
#    q=7
#    for rot in [[i,j,k] for k in range(q) for j in range(q) for i in range(q)]:
#        transrot = [13.93815383, -23.30498079, 9.44130034] + rot
#        #dimers[0].RotateBiounit(transrot[0:3], transrot[3:6])
#        print rot, dimers[0].Objective(dimers[1], transrot, biounit=True, ca=True)

