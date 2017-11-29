from collections import defaultdict
from graph import Graph
from topoltools import plot_graph
from atomdat import cov_rad
import numpy as np
from copy import copy, deepcopy

class Topology(object):

    def __init__(self, atoms, bonds=None, angles=None, dihedrals=None, impropers=None, onefour=0.5, sorted=False, ring_improper=True, bond_fudge=1.0, onlygraph=False, periodic=False):
        """
        Class for handling molecular topologies
        """

        if periodic and not all(atoms.get_pbc()):
            raise RuntimeError('Requested periodic treatment of system, but all boundary conditions set to False')

        if periodic and np.all(atoms.get_cell() == 0.):
            raise RuntimeError('Requested periodic treatment of system, but all cell dimensions set to zero')

        self.periodic = periodic
        self.atoms = atoms
        self._bonds = bonds
        self._distmat = self.calculate_distmat()
        self.onefour = onefour
        self.bond_fudge = bond_fudge

        if bonds is None:
            self._bonds = self.determine_bonds(sorted)
        else:
            self._bonds = bonds

        self.graph = Graph(vertices=range(len(atoms)), edges=self._bonds)

        if not onlygraph:
            if angles is None:
                self._angles = self.determine_angles(sorted)
            else:
                self._angles = angles

            if dihedrals is None:
                self._diheds = self.determine_dihedrals(sorted)
            else:
                self._diheds = dihedrals

            if impropers is None:
                self._improp = self.determine_impropers(sorted, ring_improper)
            else:
                self._improp = impropers

    def calculate_distmat(self):
        """
        """

        distmat = []
        alist = range(len(self.atoms))

        for i, a in enumerate(self.atoms):
            distmat.append(self.atoms.get_distances(i, alist, mic=self.periodic))

        return distmat

    def determine_bonds(self, sorted=False):
        """
        """

        bonds = []
        natoms = len(self.atoms)
        syms = self.atoms.get_chemical_symbols()

        for i in range(natoms):
            for j in range(i+1, natoms):
                if self._distmat[i][j] <= acceptable_bond_lengths(syms[i],syms[j])*self.bond_fudge:
                    l = [i, j]
                    if sorted:
                        l.sort()
                    bonds.append(l)

        return bonds

    def determine_angles(self, sorted=False):
        """
        """

        angles = []
        for vertex in self.graph:
            edges = self.graph[vertex]
            nedge = len(edges)

            for c, i in enumerate(edges):
                for j in edges[c+1:]:
                    l = [i, vertex, j]

                    if sorted:
                        l.sort()
                    angles.append(l)

        return angles

    def exclusions(self):
        gdst = self.graph.build_distance_matrix()

        natoms = len(self.atoms)
        excl = np.ones((natoms, natoms))

        for vx1 in self.graph:
            for vx2 in self.graph:

                if gdst[vx1][vx2] == 3:
                    excl[vx1][vx2] = self.onefour
                elif gdst[vx1][vx2] < 3:
                    excl[vx1][vx2] = 0.

        return excl

    def determine_dihedrals(self, sorted=False):
        """
        """

        dihed = []

        for vx1 in self.graph:
            newpath = [vx1]

            for vx2 in self.graph[vx1]:
                if vx2 not in newpath:
                    newpath = newpath[:1] + [vx2]
                else:
                    continue

                for vx3 in self.graph[vx2]:
                    if vx3 not in newpath:
                        newpath = newpath[:2] + [vx3]
                    else:
                        continue

                    for vx4 in self.graph[vx3]:
                        if vx4 not in newpath:
                            newpath = newpath[:3] + [vx4]
                        else:
                            continue

                        path = deepcopy(newpath)
                        dih = deepcopy(path)
                        dih.sort()

                        duplicate = False
                        for ndi in dihed:
                            wrk = deepcopy(ndi)
                            wrk.sort()
                            if wrk == dih:
                                #print 'duplicate'
                                duplicate = True
                                continue

                        if duplicate:
                            continue

                        if sorted:
                            path.sort()

                        dihed.append(path)

        return dihed

    def determine_impropers(self, sorted=False, ringcheck=True):
        """
        """

        improper = []
        rings = self.graph.detect_rings()
        #if at least three atoms are part of the ring then it's an improper?!?

        for vertex in self.graph:
            edges = self.graph[vertex]
            nedge = len(edges)

            #I think this is not general enough, though
            if nedge == 3:

                if ringcheck and not any([vertex in ring for ring in rings]):
                    continue

                #l = [vertex]
                l = [edges[0], edges[1], vertex, edges[2]]
                #l.extend([edges[0], edges[1], edges[2]])

                if sorted:
                    l.sort()
                improper.append(l)

        return improper

    def left_right_of_dihedral(self, dihedral):
        """
        """

        left = []
        right = []

        # a dihedral looks like this
        #          4
        #         /
        #     2--3
        #    / 
        #   1
        #
        # from all vertices, we seek the shortest path to 2 and 3
        #   if path to 2 contains 4
        #       then it's right off the dihedral
        #   elif path to 2 contains 1
        #       then it's left off the dihedral
        #   do the same with path to 3

        lchk = dihedral[0]
        rchk = dihedral[3]

        lbnd = dihedral[1]
        rbnd = dihedral[2]

        for vertex in self.graph:

            if vertex == lchk or vertex == lbnd:
                left.append(vertex)
                continue
            elif vertex == rchk or vertex == rbnd:
                right.append(vertex)
                continue

            #first check if the vertex is directly connected to dihedral:

            if vertex in self.graph[lchk] or vertex in self.graph[lbnd]:
                left.append(vertex)
                continue
            elif vertex in self.graph[rchk] or vertex in self.graph[rbnd]:
                right.append(vertex)
                continue

            ldst, lpth = self.graph.shortest_path(vertex, lbnd)
            rdst, rpth = self.graph.shortest_path(vertex, rbnd)

            #print 'left path: ', path

            if not rbnd in lpth and lbnd in rpth:
                left.append(vertex)
                continue
            elif not lbnd in rpth and rbnd in lpth:
                right.append(vertex)
                continue

            if not lbnd in rpth and not rbnd in lpth:
                #FUDO| can I do <= , = is arbitrary
                if ldst <= rdst:
                    left.append(vertex)
                    continue
                elif ldst > rdst:
                    right.append(vertex)
                    continue

            if lchk in lpth and not rchk in lpth:
                left.append(vertex)
                continue
            elif rchk in rpth and not lchk in rpth:
                right.append(vertex)
                continue

            dist, path = self.graph.shortest_path(vertex, lchk)
            #print 'left path: ', path

            if not lbnd in path and not rchk in path:
                left.append(vertex)
                continue
            else:
                dist, path = self.graph.shortest_path(vertex, rchk)
                #print 'right path: ', path

                if not rbnd in path and not lchk in path:
                    right.append(vertex)
                    continue

            raise RuntimeError("Vertex %i neither right nor left of dihedral %i %i %i %i" %(vertex, dihedral[0], dihedral[1], dihedral[2], dihedral[3]))
            #print 'not in any: ', vertex

        total = left + right
        if len(total) > len(self.atoms):
            raise RuntimeError("something went horribly in assigning left and right off of dihedral")

        return left, right

    def extract_types(self, picky=False):
        """
        """

        #naming scheme numerical for "heavy" atoms
        nametrack = defaultdict(lambda: 0)

        syms = self.atoms.get_chemical_symbols()

        fsttypes = defaultdict(list)
        fulltypes = defaultdict(list)
        shortid = {}
        loctypes = []
        shrtypes = []

        for i, a in enumerate(self.atoms):
            type = syms[i] + '#'

            if picky:
                connecting = []

                for vertex in self.graph[i]:
                    connecting.append(syms[vertex])

                    for b in self.graph[vertex]:
                        connecting[-1] += syms[b]
            else:
                #connecting = [syms[b] for b in self.graph[i]]
                #slightly more picky
                connecting = [syms[b]+str(len(self.graph[b])) for b in self.graph[i]]
                if syms[i] == 'H':
                    connecting = []

                    if len(self.graph[i]) != 1:
                        raise RuntimeError('do not know how to handle hyper-hydrogen')

                    #for vertex in self.graph[i]:
                    #    connecting.append(syms[vertex])

                    #    for b in self.graph[vertex]:
                    #        connecting[-1] += syms[b]

                    for vertex in self.graph[i]:
                        connecting.append(syms[vertex])
                        connecting[-1] += str(len(self.graph[vertex]))

            connecting.sort()

            #contract into "chemical" formula - still missing
            chem = defaultdict(lambda: 0)

            for c in connecting:
                chem[c] += 1

            l = lambda x: str(chem[x]) if chem[x] > 1 else '' 

            type += ''.join(key+l(key) for key in chem)

            if not type in loctypes: 
                shorttype = syms[i] + str(nametrack[syms[i]])
                nametrack[syms[i]] += 1
                shortid[type] = shorttype
            else:
                shorttype = shortid[type]

            shrtypes.append(shorttype)
            fsttypes[shorttype].append(i)
            fulltypes[shorttype] = type
            loctypes.append(type)

        return shrtypes

def acceptable_bond_lengths(s1, s2):
    """
    """

    return cov_rad[s1.upper()] + cov_rad[s2.upper()]
