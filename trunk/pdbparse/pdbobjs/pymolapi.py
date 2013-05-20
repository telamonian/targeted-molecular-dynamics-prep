'''
Created on Apr 1, 2013

@author: tel
'''
from tempfile import NamedTemporaryFile

from pymol import cmd

def AtomLinesToPymol(lines, name, format='pdb'):
    #tmpf = NamedTemporaryFile()
    tmpf = open('/Users/tel/Documents/workspace/md_prep/src/' + name + '_biounit.pdb', 'w')
    tmpf.write(lines)
    tmpf.flush()
    cmd.load(tmpf.name, object=name, format=format)
    tmpf.close()
    
def PymolAlign(name1, sel_string1, name2, sel_string2):
    mobile = '%s AND %s' % (name1, sel_string1)
    target = '%s AND %s' % (name2, sel_string2)
    cmd.align(mobile, target)

def PymolSave(filename, object, format):
    cmd.save(filename, selection=object, format=format)