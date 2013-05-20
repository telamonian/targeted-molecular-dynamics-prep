'''
Created on Apr 16, 2013

@author: tel
'''
import time
import numpy as np
from scipy.optimize import minimize
from collections import deque

from transmat import rotation_matrix, translation_matrix

class StruclignMixin(object):

    def PairRMSD(self, biounit=False, ca=True):
        sqdists = []
        for atom in self.IterAtoms(biounit, ca):
            sqdist = np.sum((atom.pair.xyz - atom.xyz)**2)
            if sqdist < 10000:
                sqdists.append(sqdist)
        return np.sqrt(np.mean(sqdists))
        
    def PairObjective(self, other, angvecpnt, biounit=False, ca=True):
        mat = rotation_matrix(angvecpnt[0], angvecpnt[1:4], angvecpnt[4:7])
        sqdists = []
        for atom in self.IterAtoms(biounit, ca):
            sqdists.append(np.sum((atom.pair.xyz - mat.dot(atom.xyz))**2))
        return np.sqrt(np.mean(sqdists))
    
    def PairObjectiveSplit(self, other, angvecpnt, biounit=False, ca=True):
        tmat = translation_matrix(angvecpnt[4:7])
        newcentroid = np.array(angvecpnt[4:7]) + self.Centroid(biounit, ca)
        rotmat = rotation_matrix(angvecpnt[0], angvecpnt[1:4], newcentroid.tolist())
        sqdists = []
        for atom in (atom for atom in self.IterAtoms(biounit, ca) if atom.pair!=None):
            sqdist = np.sum((atom.pair.xyz - rotmat.dot(tmat.dot(atom.xyz)))**2)
            if sqdist < 10000000:
                sqdists.append(sqdist)
        if len(sqdists) > 0:
            return np.sqrt(np.mean(sqdists))
        else:
            return 99.0

    def UnbiasedObjective(self, other, angvecpnt, biounit=False, ca=True):
        mins = []
        tmat = translation_matrix(angvecpnt[4:7])
        newcentroid = np.array(angvecpnt[4:7]) + self.Centroid(biounit, ca)
        rotmat = rotation_matrix(angvecpnt[0], angvecpnt[1:4], newcentroid.tolist())
        if self.counter % 5==0:
            oatoms = [oatom for oatom in other.IterAtoms(biounit,ca)]
            for satom in self.IterAtoms(biounit,ca):
                oatoms.sort(key=lambda atom: np.sum((atom.xyz - rotmat.dot(tmat.dot(satom.xyz)))**2))
                satom.nearby = oatoms[:10]
        self.counter += 1
        for satom in self.IterAtoms(biounit, ca):
            mins.append(np.min([np.sum((oatom.xyz - rotmat.dot(tmat.dot(satom.xyz)))**2) for oatom in satom.nearby]))
        return np.sqrt(np.mean(mins))

    def PairClose(self, other, biounit=False, ca=True):
        def Closure(angvecpnt):
            return self.PairObjective(other, angvecpnt, biounit, ca)
        return Closure
    
    def PairCloseSplit(self, other, biounit=False, ca=True):
        def Closure(angvecpnt):
            return self.PairObjectiveSplit(other, angvecpnt, biounit, ca)
        return Closure
    
    def UnbiasedClose(self, other, biounit=False, ca=True):
        def Closure(angvecpnt):
            return self.UnbiasedObjective(other, angvecpnt, biounit, ca)
        return Closure
    
    @staticmethod
    def PairCloseTrans(closure, trans):
        def Closure(rot):
            return closure(rot.tolist()+trans)
        return Closure
    
    @staticmethod
    def PairCloseRot(closure, rot):
        def Closure(trans):
            return closure(rot+trans.tolist())
        return Closure
    
    def PairMin(self, other, srange_tups, orange_tups, biounit=False, ca=True):
        initresults = []
        closure = self.PairCloseSplit(other, biounit, ca)
        ang0 = [3]
        vec0 = [.25,.75,1]
        pnt0 = ((other.Centroid(biounit, ca) - self.Centroid(biounit, ca))/2).tolist()[0:3]
        angvecpnt0 = ang0 + vec0 + pnt0
        for i in self.ThreadRange(other, srange_tups, orange_tups, biounit, ca):
            res = minimize(closure, angvecpnt0, method='TNC', bounds=((0, 6.29),(0,1),(0,1),(0,1),(-100,100),(-100,100),(-100,100)), options={'maxiter':100, 'disp':False})
            print [res.fun, res.x.tolist()]
            initresults.append([res.fun, res.x.tolist()])
#        initresults.sort()
#        results = []
#        closure = self.UnbiasedClose(other, biounit, ca)
#        for angvecpnt in [result[1] for result in initresults[:5]]:
#            self.counter = 0
#            res = minimize(closure, angvecpnt, method='TNC', bounds=((0, 6.29),(0,1),(0,1),(0,1),(-100,100),(-100,100),(-100,100)), options={'maxiter':25, 'disp':False})
#            print [res.fun, res.x.tolist()]
#            results.append([res.fun, res.x.tolist()])
        return min(initresults)
    
    def PairMinSplit(self, other, biounit=False, ca=True):
        results = []
        closure = self.PairCloseSplit(other, biounit, ca)
        for i in self.ThreadOther(other, biounit, ca):
            rot0 = [3.14,.5,.5,.5]
            trans0 = ((other.Centroid(biounit, ca) - self.Centroid(biounit, ca))/1).tolist()[0:3]
            for i in range(5):
                closetrans = self.PairCloseTrans(closure, trans0)
                res = minimize(closetrans, rot0, method='TNC', bounds=((0, 6.29),(0,1),(0,1),(0,1)), options={'maxiter':100, 'disp':False})
                rot0 = res.x.tolist()
                closerot = self.PairCloseRot(closure, rot0)
                res = minimize(closerot, trans0, method='TNC', bounds=((-100,100),(-100,100),(-100,100)), options={'maxiter':10, 'disp':False})
                trans0 = res.x.tolist()
            print [res.fun, rot0+trans0]
            results.append([res.fun, rot0+trans0])
        return min(results)
    
    def ThreadRange(self, other, srange_tups, orange_tups, biounit=False, ca=True):
        for i in range(-10,11,1):
            for atom in self.IterAtoms(biounit, ca):
                atom.pair = None
            for srange_tup, orange_tup in zip(srange_tups, orange_tups):
                for j,k in zip(range(srange_tup[1][0]+i, srange_tup[1][1]+i), range(orange_tup[1][0], orange_tup[1][1])):
                    try: 
                        self.GetChains(biounit)[srange_tup[0]].ress[j].atoms['CA'].pair = other.GetChains(biounit)[orange_tup[0]].ress[k].atoms['CA']
                    except KeyError:
                        pass
            yield 0
        
    def ThreadOther(self, other, biounit=False, ca=True):
        oatoms = deque([x for x in other.IterAtoms(biounit,ca)])
        for i in range(len(oatoms)):
            for j, (satom, oatom) in enumerate(zip(self.IterAtoms(biounit, ca), list(oatoms))):
                satom.pair = oatom
#                if j < 10:
#                    print [satom.resseq, oatom.resseq]
            yield 0
            oatoms.rotate(1)
        
    def MinRMSD(self, other, biounit=False, ca=True):
        sl, ol = len([x for x in self.IterAtoms(biounit, ca)]), len([x for x in other.IterAtoms(biounit,ca)])
        mins = np.zeros((sl,))
        oxyz = np.zeros((4,ol))
        for i, atom in enumerate(other.IterAtoms(biounit,ca)):
            oxyz[:,i] = atom.xyz[:,0]
        for i, atom in enumerate(self.IterAtoms(biounit,ca)):
            difs = oxyz - atom.xyz
            sqdists = np.array([np.sum(dif**2) for dif in difs.T])
            mins[i] = np.min(sqdists)
        return np.sqrt(np.mean(mins))
    
    def Objective(self, other, transrot, biounit=False, ca=True):
        now = time.time()
        mins = []
        mat = rotation_matrix(transrot[5], (0,0,1), self.Centroid(biounit, ca)).dot(rotation_matrix(transrot[4], (0,1,0), self.Centroid(biounit, ca))).dot(rotation_matrix(transrot[3], (1,0,0), self.Centroid(biounit, ca))).dot(translation_matrix(transrot[0:3]))
        if self.counter % 10==0:
            oatoms = [oatom for oatom in other.IterAtoms(biounit,ca)]
            for satom in self.IterAtoms(biounit,ca):
                oatoms.sort(key=lambda atom: np.sum((atom.xyz - mat.dot(satom.xyz))**2))
                satom.nearby = oatoms[:100]
#            satom_gen = self.IterAtoms(biounit,ca)
#            for i in range(20):
#                satom_gen.next()
#            test_satom = satom_gen.next()
#            print test_satom.resseq
#            for oatom in test_satom.nearby:
#                print oatom.resseq
        self.counter += 1
        for satom in self.IterAtoms(biounit, ca):
            mins.append(np.min([np.sum((oatom.xyz - mat.dot(satom.xyz))**2) for oatom in satom.nearby]))
        print 'one objective took %f' % (time.time() - now)
        print transrot
        return np.sqrt(np.mean(mins))
    
    def Close(self, other, biounit=False, ca=True):
        def Closure(transrot):
            return self.Objective(other, transrot, biounit, ca)
        return Closure
    
    @staticmethod
    def CloseTrans(closure, trans):
        def Closure(rot):
            return closure(trans+rot.tolist())
        return Closure
    
    @staticmethod
    def CloseRot(closure, rot):
        def Closure(trans):
            return closure(trans.tolist()+rot)
        return Closure
    
    def MinVec(self, other, biounit=False, ca=True):
        #bnds = ((-100, 100),(-100, 100),(-100, 100),(0,3.15),(0,3.15),(0,3.15))
        closure = self.Close(other, biounit, ca)
        trans0 = ((other.Centroid(biounit, ca) - self.Centroid(biounit, ca))/1).tolist()[0:3]
        rot0 = [4.4, 4.3, 4]
        x0 = trans0+rot0
        res = minimize(closure, x0, method='CG', bounds=((-100,100),(-100,100),(-100,100),(0,6.29),(0,6.29),(0,6.29)), options={'maxiter':10, 'disp':True, 'eps':.1, 'gtol':1})
#        for i in range(5):
#            closetrans = self.CloseTrans(closure, trans0)
#            res = minimize(closetrans, rot0,  method='TNC', bounds=((0,6.29),(0,6.29),(0,6.29)), options={'maxiter':10, 'disp':True})
#            rot0 = res.x.tolist()
#            print rot0
#            closerot = self.CloseRot(closure, rot0)
#            res = minimize(closerot, trans0,  method='Nelder-Mead', options={'maxiter':10, 'disp':True})
#            trans0 = res.x.tolist()
        print 'Final RMSD: %f' % closure(res.x)
        print res.x
        return res.x[0:3], res.x[3:6]
        #return trans0, rot0
        
    def MinVecSplit(self, other, biounit=False, ca=True):
        #bnds = ((-100, 100),(-100, 100),(-100, 100),(0,3.15),(0,3.15),(0,3.15))
        closure = self.Close(other, biounit, ca)
        trans0 = ((other.Centroid(biounit, ca) - self.Centroid(biounit, ca))/4).tolist()[0:3]
        rot0 = [4.4, 4.3, 4]
        for i in range(5):
            closetrans = self.CloseTrans(closure, trans0)
            res = minimize(closetrans, rot0,  method='TNC', bounds=((0,6.29),(0,6.29),(0,6.29)), options={'maxiter':10, 'disp':True})
            rot0 = res.x.tolist()
            print rot0
            closerot = self.CloseRot(closure, rot0)
            res = minimize(closerot, trans0,  method='Nelder-Mead', options={'maxiter':10, 'disp':True})
            trans0 = res.x.tolist()
        print 'Final RMSD: %f' % closure(res.x)
        print trans0+rot0
        return trans0, rot0
        #return trans0, rot0