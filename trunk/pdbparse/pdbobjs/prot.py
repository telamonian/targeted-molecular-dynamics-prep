'''
Created on Mar 29, 2013

@author: tel
'''
import copy, string, time
import numpy as np
from scipy.optimize import minimize
from collections import OrderedDict, deque, namedtuple

from ..data import Data
from pymolapi import AtomLinesToPymol, PymolAlign, PymolSave
from vmdapi import VMDProt
from transmat import rotation_matrix, translation_matrix
from struclign import StruclignMixin

COUNTER = 0
HETDICT = {'HOH':1,'TIP':1,'TIP3':1,'GOL':1,'RHQ':1,'KAN':1,'IMD':1,'P4P':1,'ACH':1,'TAC':1,'PUY':1,'ET':1,'M4A':1,'BER':1}

class Prot(StruclignMixin, VMDProt):
    def __init__(self, pdbparse):
        self.filename = pdbparse.pdbpath.split('/')[-1]
        self.name = self.filename.split('.')[0]
        if '.' in self.filename:
            self.format = self.filename.split('.')[1]
        else:
            self.format = 'pdb'
        self.counter = 0
        self.rules = pdbparse.rules
        self.chains = OrderedDict()
        atoms = pdbparse.rules['ATOM'].data + pdbparse.rules['HETATM'].data
        atomspec = Data(data_tuples=pdbparse.rules['ATOM'].GetDataTuples(), dtype=pdbparse.rules['ATOM'].dtype, ppm=pdbparse.rules['ATOM'].ppm, ppc=pdbparse.rules['ATOM'].ppc)
        hetspec =  Data(data_tuples=pdbparse.rules['HETATM'].GetDataTuples(), dtype=pdbparse.rules['HETATM'].dtype, ppm=pdbparse.rules['HETATM'].ppm, ppc=pdbparse.rules['HETATM'].ppc)
        for i in range(len(atoms)):
            atoms[i] = Atom(atoms[i], atomspec, hetspec)
        resi = atoms[0].resseq
        temp_atoms = []
        for atom in atoms:
            if atom.resseq != resi:
                try:
                    self.chains[temp_atoms[0].chainid].ress[resi] = Res(self, temp_atoms)
                except KeyError:
                    self.chains[temp_atoms[0].chainid] = Chain(temp_atoms[0])
                    self.chains[temp_atoms[0].chainid].ress[resi] = Res(self, temp_atoms)
                resi = atom.resseq
                temp_atoms = []
            temp_atoms.append(atom)
        #to deal with the last residue in the pdb
        try:
            self.chains[temp_atoms[0].chainid].ress[resi] = Res(self, temp_atoms)
        except KeyError:
            self.chains[temp_atoms[0].chainid] = Chain(temp_atoms[0])
            self.chains[temp_atoms[0].chainid].ress[resi] = Res(self, temp_atoms)
    
    def GetChains(self, biounit=False, target=False):
        if biounit:
            chains = self.biochains
        elif target:
            chains = self.targetchains
        else:
            chains = self.chains
        return chains.values()
    
    def IterChains(self, biounit=False, target=False):
        if biounit:
            chains = self.biochains
        elif target:
            chains = self.targetchains
        else:
            chains = self.chains
        for chain in chains.itervalues():
            yield chain
    
    def IterRess(self, biounit=False, target=False):
        for chain in self.IterChains(biounit, target):
            for res in chain.ress.itervalues():
                yield res
    
    def IterAtoms(self, biounit=False, target=False, ca=False):
        '''
        biounit=use chains from biounit, target=use chains from target, ca=use only c-alphas
        '''
        for res in self.IterRess(biounit, target):
            if ca:
                try:
                    yield res.atoms['CA']
                except KeyError:
                    try:
                        yield res.atoms["C4'"]
                    except KeyError:
                        continue
            else:
                for atom in res.atoms.itervalues():
                    yield atom
    
    def GetLength(self, biounit=False, target=False):
        lens = []
        for res in self.IterRess(biounit, target):
            lens.append(len(res.atoms))
        return np.sum(lens)
    
    def MakeBioUnit(self, biomts):
        unused_chainids = [cid for cid in string.uppercase]
        self.biochains = OrderedDict()
        for biomt in biomts:
            for chain in self.IterChains():
                if biomt.HasChain(chain.chainid):
                    tmpchain = copy.deepcopy(chain)
                    if chain.chainid in unused_chainids:
                        unused_chainids.pop(unused_chainids.index(chain.chainid))
                    else:
                        tmpchain.SetHier('chainid', unused_chainids.pop(0))
                    for res in tmpchain.ress.itervalues():
                        for atom in res.atoms.itervalues():
                            atom.xyz = biomt.mat.dot(atom.xyz)
                    self.biochains[tmpchain.chainid] = tmpchain
    
    def SepHet(self, biounit=False, newcid='Z'):
        seqresdict = {} #used to prevent residues in 'chain' Z from sharing seqres
        chains = self.chains if not biounit else self.biochains
        hetchain = Chain()
        hetchain.SetHier('chainid', newcid)
        for chain in self.IterChains(biounit):
            for key in reversed(chain.ress):
                if chain.ress[key].dtype=='HETATM':
                    tempres = chain.ress.pop(key)
                    if tempres.resseq in seqresdict:
                        for atom in tempres.IterAtoms():
                            atom.resseq = tempres.resseq*10 +1
                    else:
                        seqresdict[tempres.resseq] = 1
                    hetchain.ress['%s%s' % (key, chain.chainid)] = tempres
        chains[newcid]=hetchain
    
    def Rotate(self, angle, direction, point=None, biounit=False, target=False):
        mat = rotation_matrix(angle, direction, point)
        for atom in self.IterAtoms(biounit, target):
            atom.xyz = mat.dot(atom.xyz)
    
    @staticmethod
    def GetRotVecVec(vec1, vec2, point=None):
        rotang = np.arccos(np.dot(vec1, vec2)/(np.sqrt(np.dot(vec1, vec1))*np.sqrt(np.dot(vec2, vec2))))
        rotax = np.cross(vec1, vec2)
        return rotation_matrix(rotang, rotax, point)
    
    def RotateSplit(self, angle, direction, point=None, biounit=False, ca=True):
        if point==None:
            point=[0,0,0]
        tmat = translation_matrix(point)
        newcentroid = np.array(point) + self.Centroid(biounit, ca)
        rotmat = rotation_matrix(angle, direction, newcentroid.tolist())
        for atom in self.IterAtoms(biounit):
            atom.xyz = rotmat.dot(tmat.dot(atom.xyz))
    
    def RotateBiounit(self, trans, rot):
        mat = rotation_matrix(rot[2], (0,0,1), self.Centroid(biounit, ca)).dot(rotation_matrix(rot[1], (0,1,0), self.Centroid(biounit, ca))).dot(rotation_matrix(rot[0], (1,0,0), self.Centroid(biounit, ca))).dot(translation_matrix(trans[0:3]))
        for atom in self.IterAtoms(biounit=True):
            atom.xyz = mat.dot(atom.xyz)
    
    def Centroid(self, biounit=False, target=False, ca=False):
        x,y,z = [],[],[]
        for atom in self.IterAtoms(biounit, ca):
            x.append(atom.x)
            y.append(atom.y)
            z.append(atom.z)
        return np.array((np.mean(x), np.mean(y), np.mean(z)))
    
    def GetOrientTransmat(self, biounit=False, target=False):
        cen = np.zeros((3,1))
        cenlist = self.Centroid(biounit, target)
        cen[0:3] = [[cenlist[0]],[cenlist[1]],[cenlist[2]]]
        tmat = translation_matrix(-self.Centroid(biounit, target))
        max = 0
        farthestxyz = None
        for atom in self.IterAtoms(biounit, target):
            dist = np.sum((atom.xyz[0:3] - cen)**2)
            if dist > max:
                farthestxyz = atom.xyz[0:3] - cen
                max = dist
        firstrotax = np.cross(farthestxyz.transpose().tolist()[0], np.array([1,0,0]))
        firstrotang = np.arccos(np.dot(farthestxyz.transpose().tolist()[0], np.array([1,0,0]))/(np.sqrt(np.dot(farthestxyz.transpose().tolist()[0], farthestxyz.transpose().tolist()[0]))))
        firstrotmat = rotation_matrix(firstrotang, firstrotax, self.Centroid(biounit, target))
        max = 0
        firsttransmat = tmat.dot(firstrotmat)
        for atom in self.IterAtoms(biounit, target):
            dist = sum(firsttransmat.dot(atom.xyz)[1:3]**2)
            if dist > max:
                secfarthestyz = firsttransmat.dot(atom.xyz)[0:3]
                max = dist
        secfarthestyz[0] = 0
        secondrotax = np.array([1,0,0])
        secondrotang = np.arccos(np.dot(secfarthestyz.transpose().tolist()[0], np.array([0,1,0]))/np.sqrt(np.dot(secfarthestyz.transpose().tolist()[0],secfarthestyz.transpose().tolist()[0])))
        secondrotmat = rotation_matrix(secondrotang, secondrotax, self.Centroid(biounit, target))
        return secondrotmat.dot(firsttransmat)
    
    def Center(self, biounit=False, target=False):
        tmat = translation_matrix(-self.Centroid(biounit, target))
        self.Transform(tmat, biounit, target)
    
    def Orient(self, others=None, biounit=False, target=False):
        pca = self.PCA(biounit, target)
        rotmat = self.GetRotVecVec(pca[0][2], np.array((1,0,0)), self.Centroid(biounit, target))
        self.Transform(rotmat, biounit, target)
        rotmat = self.GetRotVecVec(pca[0][0], np.array((0,1,0)), self.Centroid(biounit, target))
        self.Transform(rotmat, biounit, target)
        if others:
            for other in others:
                other.Transform(rotmat, biounit, target)
    
    def PCA(self, biounit=False, target=False):
        '''
        taken from http://glowingpython.blogspot.it/
        '''
        numpc = 3
        #creating array of atom coordinates
        coords = np.zeros((self.GetLength(biounit, target),3))
        for i, atom in enumerate(self.IterAtoms(biounit, target)):
            coords[i,:] = atom.xyz[0:3].transpose()
        # computing eigenvalues and eigenvectors of covariance matrix
        M = (coords - np.mean(coords.T,axis=1)).T # subtract the mean (along columns)
        [latent,coeff] = np.linalg.eig(np.cov(M))
        p = np.size(coeff,axis=1)
        idx = np.argsort(latent) # sorting the eigenvalues
        idx = idx[::-1]       # in ascending order
        # sorting eigenvectors according to the sorted eigenvalues
        coeff = coeff[:,idx]
        latent = latent[idx] # sorting eigenvalues
        if numpc < p or numpc >= 0:
            coeff = coeff[:,range(numpc)] # cutting some PCs
        score = np.dot(coeff.T,M) # projection of the data in the new space
        return coeff,score,latent 
    
    def Transform(self, transmat, biounit=False, target=False):
        for atom in self.IterAtoms(biounit, target):
            atom.xyz = transmat.dot(atom.xyz)
    
    def GetTarget(self, other, biounit):
        for satom in self.IterAtoms(biounit, ca=True):
            mini = 999.0
            for chain in other.targetchains.itervalues():
                for oatom in chain.IterAtoms(ca=True):
                    dist = satom.Dist(oatom)
                    if dist < mini and (oatom.resseq - 10 < satom.resseq < oatom.resseq + 10):
                        satom.target = oatom
                        mini = dist
    
    def MakeTargetPDB(self, other, outpath, biounit=False, ca=False):
        self.targetchains = OrderedDict()
        for chain in self.IterChains(biounit):
            tmpchain = copy.deepcopy(chain)
            for atom in tmpchain.IterAtoms():
                atom.occupancy = 0.
            self.targetchains[chain.chainid] = tmpchain
        other.GetTarget(self, biounit)
        for oatom in other.IterAtoms(biounit, ca=True):
            oatom.target.xyz = oatom.xyz
            oatom.target.occupancy = 1.
        self.Save(outpath, target=True, ca=ca)
    
    def RemoveNamelist(self, namelist, biounit=False):
        for chain in self.IterChains(biounit):
            chain.RemoveNamelist(namelist)
            
    def RemoveFarWater(self, cutoff=2.1, biounit=False, target=False):
        for chain in self.IterChains(biounit):
            chain.RemoveFarWater(self.GetChains(biounit, target), cutoff)
    
    def RenumberNegs(self, biounit=False, target=False):
        for chain in self.IterChains(biounit, target):
            foundneg = False
            for res in chain.IterRess():
                if res.resseq < 1:
                    foundneg = True
                    break
            if foundneg:
                counter = 1
                for res in chain.IterRess():
                    res.resseq = counter
                    for atom in res.IterAtoms():
                        atom.resseq = counter
                    counter +=1
    
    def GetDims(self, biounit=False, target=False):
        min = np.zeros((1,3))
        min[0,:] = 9999999
        max = -min
        for atom in self.IterAtoms(biounit, target):
            for i, dim in enumerate(atom.xyz[:3]):
                if dim > max[0,i]:
                    max[0,i] = dim
                elif dim < min[0,i]:
                    min[0,i] = dim
        return min, max, max-min, (max+min)/2.0
        
    
    def Print(self, biounit=False, target=False, ca=False, posocc=False):
        out = ''
        if posocc:
            for atom in self.IterAtoms(biounit, target, ca):
                if atom.occupancy > 0:
                    out += '%s\n' % atom.Print()#self.rules['ATOM'], hetspec=self.rules['HETATM'])
        else:
            for atom in self.IterAtoms(biounit, target, ca):
                out += '%s\n' % atom.Print()#atomspec=self.rules['ATOM'], hetspec=self.rules['HETATM'])
        return out
    
    def Save(self, fpath, biounit=False, target=False, ca=False, posocc=False):
        '''
        fpath=filepath, biounit=use chains from biounit, target=use chains from target, ca=use only c-alphas, posocc=use only atoms with positive occupancy
        '''
        with open(fpath, 'w') as f:
            for line in self.Print(biounit, target, ca, posocc):
                f.write(line)
    
    def SaveChains(self, fpath, biounit=False, target=False, ca=False, posocc=False):
        for chain in self.IterChains(biounit, target):
            with open(fpath.split('.')[0]+('_%s' % chain.chainid)+'.pdb', 'w') as f:
                for line in chain.Print(ca, posocc):
                    f.write(line)
    
    def PymolLoad(self, biounit=False):
        AtomLinesToPymol(self.Print(biounit), name=self.name, format=self.format)
        
    def PymolAlign(self, other, sel1, sel2):
        PymolAlign(name1=self.name, sel_string1=sel1, name2=other.name, sel_string2=sel2)
    
    def PymolSave(self, filename=None):
        if filename == None:
            filename = self.filename
        PymolSave(filename, self.name, self.format)
        
    @staticmethod
    def MergeSolvTarget(tpath, spath):
        '''
        combines non-waters from the _target pdb file with waters from the
        _vmd_wat_ion pdb file
        '''
        with open(tpath) as tf:
            tlines = tf.readlines()
        with open(spath) as sf:
            slines = sf.readlines()
        outpath = ''.join(tpath.split('.')[:-1] + ['_wat_ion.'] + [tpath.split('.')[-1]])
        with open(outpath, 'w') as out:
            for i, line in enumerate(tlines):
                '''this loop writes non-waters from the _target pdb'''
                if line[17:20] == 'TIP':
                    break
                out.write(line)
            for line in slines[i+1:]:
                '''this loop writes waters from the _vmd_wat_ion pdb'''
                line = line[:54] + '  0.00' + line[60:] #makes sure the water occupancy is set to 0.00
                out.write(line)
                
class Chain(object):
    def __init__(self, atom=None):
        if atom!=None:
            self.chainid = atom.chainid
        self.ress = OrderedDict()
        
    def IterRess(self):
        for res in self.ress.itervalues():
            yield res
    
    def IterAtoms(self, ca=False):
        for res in self.IterRess():
            if ca:
                try:
                    yield res.atoms['CA']
                except KeyError:
                    try:
                        yield res.atoms["C4'"]
                    except KeyError:
                        continue
            else:
                for atom in res.atoms.itervalues():
                    yield atom
    
    def RemoveNamelist(self, namelist):
        for key in reversed(self.ress):
            if self.ress[key].resname in namelist:
                del self.ress[key]
    
    def RemoveFarWater(self, chains, cutoff=2.1):
        for key in reversed(self.ress):
            if self.ress[key].resname=='HOH':
                min = 999
                preserve = False
                for chain in chains:
                    for atom in (atom for atom in chain.IterAtoms() if atom.dtype=='ATOM'):
                        if atom.Dist(self.ress[key].atoms['O']) < min:
                            min = atom.Dist(self.ress[key].atoms['O'])
                        if atom.Dist(self.ress[key].atoms['O']) <= cutoff:
                            preserve = True
                            break
                    if preserve:
                        break
                if not preserve:
                    del self.ress[key]
                print min
                    
    def SetHier(self, name, val):
        self.__setattr__(name, val)
        for res in self.ress.itervalues():
            res.SetHier(name, val)
            
    def Print(self, ca=False, posocc=False):
        out = ''
        if posocc:
            for atom in self.IterAtoms(ca):
                if atom.occupancy > 0:
                    out += '%s\n' % atom.Print()#self.rules['ATOM'], hetspec=self.rules['HETATM'])
        else:
            for atom in self.IterAtoms(ca):
                out += '%s\n' % atom.Print()#atomspec=self.rules['ATOM'], hetspec=self.rules['HETATM'])
        return out

class Res(object):
    def __init__(self, prot, atoms):
        self.prot = prot
        self.dtype = atoms[0].dtype
        self.resseq = int(atoms[0].resseq)
        self.resname = atoms[0].resname
        self.chainid = atoms[0].chainid
        self.atoms = OrderedDict()
        for atom in atoms:
            self.atoms[atom.name] = atom
            
    def IterAtoms(self, ca=False):
        if ca:
            try:
                yield self.atoms['CA']
            except KeyError:
                try:
                    yield self.atoms["C4'"]
                except KeyError:
                    pass
        else:
            for atom in self.atoms.itervalues():
                yield atom
    
    def GetNormVec(self):
        if self.resname=='DA' or self.resname=='DG':
            vec1 = self.atoms['N9'].GetVecTo(self.atoms['C8'])
            vec2 = self.atoms['N9'].GetVecTo(self.atoms['C4'])
        elif self.resname=='DT' or self.resname=='DC':
            vec1 = self.atoms['N1'].GetVecTo(self.atoms['C6'])
            vec2 = self.atoms['N1'].GetVecTo(self.atoms['C2'])
        else:
            raise
        return np.cross(vec1, vec2)
    
    def GetRoll(self, other):
        norm1 = self.GetNormVec()
        norm2 = other.GetNormVec()
        return np.arccos(np.dot(norm1,norm2)/(np.sqrt(np.dot(norm1, norm1))*np.sqrt(np.dot(norm2, norm2))))
    
    def GetNextRes(self):
        return self.prot.chains[self.chainid].ress[self.resseq+1]
    
    def SetHier(self, name, val):
        self.__setattr__(name, val)
        for atom in self.atoms.itervalues():
            atom.SetHier(name, val)
    
    def Build(self, fpath='agluc.txt'):
        Ic = namedtuple('Ic', ('i','j','k','l','r1','t1','p','t2','r2'))
        ics = []
        imps = []
        f = open(fpath)
        for line in f:
            line = line.split()
            if '*' in line[3]:
                imps.append(Ic([token.strip('*') for token in line[1:]]))
            else:
                ics.append(Ic(line[1:]))
        f.close()
        first = ics[0]
        for index in 'ijk':
            pass
        for ic in ics:
            pass
    
    @staticmethod
    def Header():
        return '%s\t%s\t%s' % ('base_number', 'base_type', 'roll_value (radians)')
    
    def __str__(self):
        return '%s\t%s\t%.3f' % (self.resseq, self.resname, self.GetRoll(self.GetNextRes()))
        
class Atom(object):
    def __init__(self, atom, atomspec, hetspec):
        self.atomspec = atomspec
        self.hetspec = hetspec
        self.xyz = np.zeros((4,1))    #assemble xyz vector into this array
        self.xyz[3,0] = 1    #use 4x1 matrix for coordinates to allow for use of translate+rotate matrices
        if atom.resName in HETDICT:
            self.dtype = 'HETATM'
        else:
            self.dtype = atom.dtype
        for tup in atom.GetDataTuples():
            i = 'xyz'.find(tup[0])
            if i >= 0:
                self.xyz[i,0] = float(tup[1])
            else:
                try:    #assume that datatype is [int, float, str], falling to the next case if an error occurs
                    self.__setattr__(tup[0].lower(), int(tup[1]))
                except ValueError:
                    try:
                        self.__setattr__(tup[0].lower(), float(tup[1]))
                    except ValueError:
                        self.__setattr__(tup[0].lower(), tup[1])

    def GetVecTo(self, other):
        return np.array([other.x - self.x, other.y - self.y, other.z - self.z])
    
    def Dist(self, other):
        return np.sqrt(np.sum(self.GetVecTo(other)**2))
    
    def SetHier(self, name, val):
        self.__setattr__(name, val)
    
    def Print(self):
        '''
        hacky (buggy?) printer of ATOM line that matches PDB spec
        doesn't respect proper float format
        '''
        spec = self.atomspec if self.dtype=='ATOM' else self.hetspec
        out = '%-*s' % (spec.ppc[1]-spec.ppc[0], spec.ppm)
        for name in spec.data_list:
            try:
                pad = spec.__getattribute__(name)[0] - last_end
                last_end = spec.__getattribute__(name)[1]
            except NameError:
                pad = 0
                last_end = spec.__getattribute__(name)[1]
            if name=='name':  #this conditional to deal with the stupid bastard cases that PDB seems to include for the hell of it
                align = '-'
                if len(self.__getattribute__(name)) < 4:
                    extra = ' '
            elif name in 'xyz':
                kind = '.3f'
            elif name=='occupancy' or name=='tempfactor':
                kind = '.2f'
            else:
                extra = ''
                kind = 's'
                align = ''
            format = '%%*s%s%%%s*%s' % (extra, align, kind)
            out += format % (pad, '', spec.__getattribute__(name)[1] - spec.__getattribute__(name)[0] - len(extra), self.__getattribute__(name.lower()))
        return out
    
    @property
    def x(self):
        return self.xyz[0,0]
    @property
    def y(self):
        return self.xyz[1,0]
    @property
    def z(self):
        return self.xyz[2,0]

