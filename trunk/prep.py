'''
Created on Mar 17, 2013

@author: tel
'''

from pdbparse.transmat import rotation_matrix

class axis(object):
    def __init__(self, abc=(0,0,0), xyz=(0,0,0)):
        self.abc = abc
        self.a = abc[0]
        self.b = abc[1]
        self.c = abc[2]
        self.xyz = xyz 
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]

    def axes_align_mat(self, other):
        