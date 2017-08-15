#!/usr/bin/env python

from moltop.topol import Topology
from ase.io import read, write
from ase import Atoms
import numpy as np

atoms = read('c4c1im_pack.xyz')
top = Topology(atoms, sorted=True)
#print top._bonds
#print top._angles

f = open('bonds.txt')

blist = []

for line in f.readlines():
    nds = map(int, line.split()[:2])
    nds[0] -= 1
    nds[1] -= 1
    nds.sort()
    blist.append(nds)

f.close()
    
alist = []
f = open('angles.txt')
for line in f.readlines():
    nds = map(int, line.split()[:3])
    nds[0] -= 1
    nds[1] -= 1
    nds[2] -= 1
    nds.sort()
    alist.append(nds)
    
f.close()
    
print('Checking bonds')
fucked = False
for b in top._bonds:
    if b not in blist:
        fucked = True

if len(top._bonds) != len(blist):
    fucked = True

if fucked:
    print('Bonds are fucked')
else:
    print('Bonds are fine')

fucked = False
        
print('Checking angles')
for a in top._angles:
    if a not in alist:
        fucked = True

if len(top._angles) != len(alist):
    fucked = True

if fucked:
    print('Angles are fucked')
else:
    print('Angles are fine')

f = open('dihedrals.txt')

dlist = []
for line in f.readlines():
    nds = map(int, line.split()[:4])
    nds[0] -= 1
    nds[1] -= 1
    nds[2] -= 1
    nds[3] -= 1
    nds.sort()
    dlist.append(nds)
    
fucked = False

dall = top._diheds
dall.extend(top._improp)

for d in dlist:
    if d not in dall:
        print d
        fucked = True

if len(dlist) != len(dall):
    fucked = True

if fucked:
    print('Dihedrals are fucked')
else:
    print('Dihedrals are fine')

rings = top.graph.detect_rings()
print 'rings are: ', rings

top = Topology(atoms, sorted=False)

ndih = 20
ndih = 0
ndih = 10
ndih = 49
ndih = 39

dihed = top._diheds[ndih]
print 'my dihedral is: ', dihed
left, right = top.left_right_of_dihedral(dihed)

#write('right.xyz', top.atoms[right])
#write('left.xyz', top.atoms[left])
#write('dihed.xyz', top.atoms[top._diheds[ndih]])

strt = top.atoms.get_dihedral(dihed)

traj = []

ltrj = []
rtrj = []
dtrj = []

for ngl in np.linspace(strt, strt + 2.*np.pi, 11):
    cpy = top.atoms.copy()
    cpy.set_dihedral(top._diheds[ndih], ngl, indices=right)
    #cpy.set_dihedral(top._diheds[ndih], ngl, indices=left)
    traj.append(cpy)
    ltrj.append(cpy[left])
    rtrj.append(cpy[right])
    dtrj.append(cpy[dihed])

write('torsion.xyz', traj)

write('right.xyz', rtrj)
write('left.xyz', ltrj)
write('dihed.xyz', dtrj)
