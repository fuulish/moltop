from collections import defaultdict
from graph import Graph
from topoltools import plot_graph
from atomdat import cov_rad
import matplotlib.pyplot as plt
import numpy as np
from copy import copy, deepcopy

class Topology(object):

    def __init__(self, atoms, bonds=None, angles=None, dihedrals=None, impropers=None, onefour=0.5, sorted=False, ring_improper=False, bond_fudge=1.0, onlygraph=False):
        """
        Class for handling molecular topologies
        """

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

        #plot_graph(atoms.positions, g)
        #plt.show()

    def calculate_distmat(self):
        """
        """

        distmat = []
        alist = range(len(self.atoms))

        for i, a in enumerate(self.atoms):
            distmat.append(self.atoms.get_distances(i, alist))

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
            #if nedge == 2:
            #    l = [edges[0], vertex, edges[1]]
            #    #l.extend(edges)
            #    if sorted:
            #        l.sort()
            #    angles.append(l)
            #elif nedge > 2:
            #print 'new: ', edges
            #if nedge >= 2:
            for c, i in enumerate(edges):
                for j in edges[c+1:]:
                    l = [i, vertex, j]
                    #print 'actual angle: ', l
                    #l.extend([i,j])
                    if sorted:
                        l.sort()
                    angles.append(l)

        #print angles

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
        #path = self.graph.build_path_matrix()
        #rings = self.graph.find_rings()

        #again, problems with ring structures,
        #because there might be a dihedral at a longer path than the shortest one that will hence not be recognized
        #i.e., determine all paths, and if more than one path, check that one as well.
        #create second graph, that removes vertices contained in currently shortest path and 
        #   check connectivity (if disconnected, only one path)
        #   otherwise keep calculating shortest path, and recalculate new shortest path
        # if ring structure -> two paths, if more than two paths, then something is weird

        #for vx1 in self.graph:
        #    one = [vx1]
        #    for vx2 in self.graph[vx1]:
        #        if vx2 in one:
        #            continue
        #        one.append(vx2)
        #        for vx3 in self.graph[vx2]:
        #            if vx3 in one:
        #                continue
        #            one.append(vx3)
        #            for vx4 in self.graph[vx3]:
        #                if vx4 in one:
        #                    continue
        #                one.append(vx4)

        #                dih = deepcopy(one)
        #                dih.sort()

        #                for ndi in dihed:
        #                    wrk = deepcopy(ndi)
        #                    wrk.sort()
        #                    if wrk == dih:
        #                        #print 'duplicate'
        #                        duplicate = True
        #                        continue

        #                if sorted:
        #                    one.sort()

        #                dihed.append(one)

        for vx1 in self.graph:
            for vx2 in self.graph:
                if vx1 == vx2:
                    continue

                paths = self.graph.find_all_paths(vx1, vx2)

                #print paths
                for path in paths:
                    if len(path) == 4:
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

    def determine_impropers(self, sorted=False, ringcheck=False):
        """
        """

        improper = []
        rings = self.graph.detect_rings()
        #if at least three atoms are part of the ring then it's an improper?!?

        for vertex in self.graph:
            edges = self.graph[vertex]
            nedge = len(edges)

            if nedge == 3:
                l = [vertex]
                l.extend([edges[0], edges[1], edges[2]])

                #ringcheck, only do impropers if in ring structure?
                if ringcheck:
                    cnt = []
                    for ring in rings:
                        cnt.append(0)
                        for e in l:
                            if e in ring:
                                cnt[-1] += 1

                    if all([c < 3 for c in cnt]):
                        continue

                if sorted:
                    l.sort()
                improper.append(l)

            ##as necessary condition need all edges at least have 2 edges as well
            ##probably not sufficient though
            ##should check if it's actually a ring structure

            #if nedge > 2:
            ##if nedge == 3:
            #    for i, e1 in enumerate(edges):
            #        for j, e2 in enumerate(edges[i+1:]):
            #            for k, e3 in enumerate(edges[j+i+2:]):
            #                print ''
            #                print vertex+1, e1+1, e2+1, e3+1
            #                print self.graph[e1]
            #                print self.graph[e2]
            #                print self.graph[e3]
            #                nbtruth = [len(self.graph[e1]) > 1, len(self.graph[e2]) > 1, len(self.graph[e3]) > 1]

            #                
            #                conn = 0
            #                for mxcon in nbtruth:
            #                    if mxcon:
            #                        conn += 1

            #                if conn > 1:
            #                    continue

            #                #print nbtruth, all(nbtruth)
            #                #if not all(nbtruth):
            #                #    continue

            #                l = [vertex, e1, e2, e3]

            #                if sorted:
            #                    l.sort()

            #                improper.append(l)


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

def acceptable_bond_lengths(s1, s2):
    """
    """

    return cov_rad[s1] + cov_rad[s2]
