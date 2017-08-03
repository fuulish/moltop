#!/usr/bin/env python

from moltop.topol import Topology
from ase.io import read, write
from ase import Atoms
import numpy as np

atoms = read('gdm.xyz')
top = Topology(atoms, sorted=False)

ndih = 0

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
    cpy.set_dihedral(top._diheds[ndih], ngl, indices=left)
    traj.append(cpy)
    ltrj.append(cpy[left])
    rtrj.append(cpy[right])
    dtrj.append(cpy[dihed])

write('torsion.xyz', traj)

write('right.xyz', rtrj)
write('left.xyz', ltrj)
write('dihed.xyz', dtrj)
