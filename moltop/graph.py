# IPython log file

from collections import defaultdict
from copy import copy, deepcopy

#FUDO| check everywhere whether to deepcopy...

class Graph(object):
    def __init__(self, vertices=None, edges=None, dists=None, check_connectivity=True):

        self.verts = set()

        #sets are nice for this because no duplicates, but they are unordered
        #could use collections.orderedict if we really need the O(1) access
        self.edges = defaultdict(list)
        self.dists = defaultdict(list)

        #self.paths = defaultdict(list)
        #self.pdist = defaultdict(list)

        for vertex in vertices:
            self.add_vertex(vertex)

        if dists is None:
            dists = len(edges)*[1,]

        for edge, dist in zip(edges, dists):
            self.add_edge(*edge, dist=dist)

        #if dists is None:
        #    for vertex in self.verts:
        #        self.dists[vertex] = defaultdict(list)
        #        for edge in self[vertex]:
        #            dists[vertex][edge].append(1)

        if check_connectivity:
            self.check_connectivity()

    def check_connectivity(self):

        if not self.connected():
            raise RuntimeError("There are disconnected vertices")

    def connected(self):

        nedge = 0

        for vertex in self.verts:
            if len(self.edges[vertex]) == 1:
                nedge += 1
            elif len(self.edges[vertex]) == 0:
            #elif not vertex in self.edges:
                return False

#FUX| guessing, double-check
        #if nedge > len(self) - 2:
        #    return False
        #    #raise RuntimeError("There are disconnected vertices")

        # at least exclude the ones already check, only check if from one end can get to all other ends
        for i, vx1 in enumerate(self.verts):
            if i > 0:
                break
            for vx2 in self.verts:
                if vx1 == vx2:
                    continue

                dist, path = self.shortest_path(vx1, vx2)

                if dist == -1:
                    return False

        return True

    def add_vertex(self, vertex):
        #don't need to check double-counting, because verts is set
        self.verts.add(vertex)
    
    def add_edge(self, from_vert, to_vert, dist=1):

        #if not to_vert in self.edges[from_vert] and not from_vert in self.edges[to_vert]:
        #    self.edges[from_vert].append(to_vert)
        #    self.dists[from_vert].append(dist)

        #    self.edges[to_vert].append(from_vert)
        #    self.dists[to_vert].append(dist)
        #else:
        #    raise RuntimeError('Edge already exists')

        if not to_vert in self.edges[from_vert]:
            self.edges[from_vert].append(to_vert)
            self.dists[from_vert].append(dist)
        else:
            raise RuntimeError('Edge already exists')

        if not from_vert in self.edges[to_vert]:
            self.edges[to_vert].append(from_vert)
            self.dists[to_vert].append(dist)
        else:
            raise RuntimeError('Edge already exists')

    def remove_vertex(self, remove, discard_disconnected=False, split_graph=False, check_connectivity=True):
        """
        """

        #FUDO| create copy of self remove vertex and associated edges
        #FUDO| check connectivity and raise corresponding error if graph disconnected

        graph = deepcopy(self)

        #remove the vertex
        if remove in graph.verts:
            graph.verts.remove(remove)

        #remove associated edges that contain remove
        for vertex in graph:
            if remove in graph[vertex]:
                graph[vertex].remove(remove)

        #if RuntimeError is raised, nothing happened to self

        if discard_disconnected:
            graph.remove_disconnected_vertices()

        #if not graph.connected() and split_graph:
        #FUX| I don't really care if they're disconnected or not, should I?
        #FUX| if it's one graph, it'll return a list of one graph
        if split_graph:
            return graph.split_graph()
        else:
            if check_connectivity:
                graph.check_connectivity()

        #FUDO| or do we want to return the result?

        return graph
        #self = deepcopy(graph)

    def split_graph(self):

        graphs = []
        vertices = None
        edges = None
        dists = None

        #print 'new split graph'

        for vx1 in self.verts:

            if any([vx1 in graph.verts for graph in graphs]):
                #print 'already contained'
                continue

            vertices = [vx1]
            edges = [[vx1, e] for e in self[vx1]]
            dists = [d for d in self.dists[vx1]]

            for vx2 in self.verts:
                if vx1 == vx2:
                    continue

                dist, path = self.shortest_path(vx1, vx2)

                if dist == -1:
                    continue
                else:
                    vertices.append(vx2)
                    new_edges = [sorted([vx2, e]) for e in self[vx2]]
                    new_dists = [d for d in self.dists[vx2]]

                    #avoid duplicates, ordering doesn't matter because edge is added to both vertices
                    for dst, edg in zip(new_dists, new_edges):
                        if not edg in edges:
                            edges.append(edg)
                            dists.append(dst)


            #print 'VERTICES: ', vertices
            #print 'EDGES: ', edges
            g = Graph(vertices=vertices, edges=edges, dists=dists, check_connectivity=False)
            graphs.append(g)

        #for g in graphs:
        #    print g.verts

        #print 'in total ', len(graphs), ' new graphs'
        return graphs

    def remove_edge(self, from_vertex, to_vertex, discard_disconnected=False, check_connectivity=True):

        graph = deepcopy(self)

        graph.edges[from_vertex].remove(to_vertex)
        graph.edges[to_vertex].remove(from_vertex)

        if discard_disconnected:
            graph.remove_disconnected_vertices()

        if check_connectivity:
            graph.check_connectivity()

        return graph

    def remove_disconnected_vertices(self):
        graph = deepcopy(self)

        for vertex in graph:
            #print vertex, self[vertex]
            if len(graph[vertex]) == 0:
                self.verts.remove(vertex)

    def detect_rings(self):
        """
        """

        rings = []

        for vx1 in self:
            edges = self[vx1]
            nedge = len(edges)

            if nedge < 2:
                continue

            for vx2 in self:

                #exclude if the same
                if vx1 == vx2:
                    continue

                paths = self.find_all_paths(vx1, vx2)

                rng = []

                if len(paths) > 1:
                    for p in paths:
                        rng.extend(p)

                    rng = set(rng)

                    if rng in rings:
                        continue
                    else:
                        rings.append(rng)

        return rings

    def find_all_paths(self, from_vertex, to_vertex):
        """
        """

        if from_vertex == to_vertex:
            raise RuntimeError('Finding paths between the same vertices does not work right now')

        graph = deepcopy(self)
        allpaths = []

        while len(graph):

            #I don't think we need this
            #if not from_vertex in graph.verts or not to_vertex in graph.verts:
            #    break

            dist, path = graph.shortest_path(from_vertex, to_vertex)

            allpaths.append(path)

            #if it's an edge
            if len(path) == 2:
                graph = graph.remove_edge(path[0], path[1], discard_disconnected=True, check_connectivity=False)
            else:
                for vertex in path[1:-1]:
                    graph = graph.remove_vertex(vertex, discard_disconnected=True, check_connectivity=False)

            if not from_vertex in graph.verts or not to_vertex in graph.verts:
                break

            dist1, path = graph.shortest_path(from_vertex, to_vertex)
            dist2, path = graph.shortest_path(to_vertex, from_vertex)

            if all([d < 0 for d in [dist1, dist2]]):
                break
            else:
                graphs = graph.split_graph()

                for g in graphs:
                    if to_vertex in g.verts and from_vertex in g.verts:
                        graph = g

        return allpaths

    def __getitem__(self, key):
        if not key in self.edges:
            raise IndexError("Vertex not in graph")

        return self.edges[key]#, self.dists[key]

    def __len__(self):
        return len(self.verts)

    def __add__(self, graph):

        newgraph = copy(self)
        #newgraph = deepcopy(self)

        #first, add all vertices
        for vertex in graph.verts:
            #print vertex
            if not vertex in newgraph.verts:
                newgraph.add_vertex(vertex)

        #second, copy all edges
        for vertex in graph.verts:
            #because we don't know if all vertices have edges? (they should have at least on edges, otherwise disconnected
            #this shouldn't happen, because we check connectivity every time we initialize graph
            if not vertex in graph.edges:
                continue

            for edge in graph.edges[vertex]:
                if not edge in newgraph.edges[vertex]:
                    newgraph.add_edge(vertex, edge)

        return newgraph

    def __iter__(self):
        return iter(self.verts)

    def __eq__(self, other):
        if not self.verts == other.verts:
            return False
        else:
            for vertex in self:
                l1 = self[vertex]
                l2 = other[vertex]

                l1.sort()
                l2.sort()

                if not l1 == l2:
                    return False
                
                l1 = self.dists[vertex]
                l2 = other.dists[vertex]

                l1.sort()
                l2.sort()

                if not l1 == l2:
                    return False

        return True

    def __ne__(self, other):
        """
        """

        iseq = self.__eq__(other)

        return not iseq

    def distance(self, from_vert, to_vert):
        #if all distances were one:

        #if self.paths[from_vert] is None:
        dist, prvs = self.shortest_path(from_vert, to_vert)

        #for prev in prvs:
        #    path = spell_out_path(prvs, from_vert, to_vert)

        return dist[to_vert]

    def build_path_matrix(self):
        """
        """

        pathmat = defaultdict(list)

        for vx1 in self:
            dists, prvs = self.shortest_path(vx1)
            pathmat[vx1] = defaultdict(list)
            for vx2 in self:
                pathmat[vx1][vx2] = self.spell_out_path(prvs, vx1, vx2)

        #self.distmat = copy(distmat)
        return pathmat

    def build_distance_matrix(self):
        """
        """

        distmat = defaultdict(list)

        for vx1 in self:
            dists, prvs = self.shortest_path(vx1)
            distmat[vx1] = defaultdict() #could just use {}
            for vx2 in self:
                distmat[vx1][vx2] = dists[vx2]

        #self.distmat = copy(distmat)
        return distmat

    def edge_length(self, from_vert, to_vert):
        """
        """

        edges = self.edges[from_vert]
        ndx = edges.index(to_vert)

        return self.dists[from_vert][ndx]

    def shortest_path(self, from_vert, to_vert=None):
        return Dijkstra(self, from_vert, to_vert)

    @staticmethod
    def spell_out_path(prev, from_vert, to_vert):
        return Dijkstra_path(prev, from_vert, to_vert)

def Dijkstra(graph, source, destination=None):
    """
    """

    unvisited = copy(graph.verts)

    dist = {}
    prev = {}

    inf = float('inf')

    for vertex in unvisited:
        dist[vertex] = inf
        prev[vertex] = -1

    dist[source] = 0

    while len(unvisited) != 0:
        #u = dist.index(min(dist))
        #unvisited.remove(u)

        mindist = inf
        minvert = -1
        for vertex in unvisited:
            if dist[vertex] < mindist:
                mindist = dist[vertex]
                minvert = vertex

        if minvert == -1:
            return -1, []

        unvisited.remove(minvert)

        if destination is not None:
            if minvert == destination:
                break

        for neighbor in graph[minvert]:
            alt = mindist + graph.edge_length(minvert, neighbor)
            if alt < dist[neighbor]:
                dist[neighbor] = alt
                prev[neighbor] = minvert

    if destination is not None:
        path = Dijkstra_path(prev, source, destination)
        return dist[destination], path
    else:
        return dist, prev

def Dijkstra_path(prev, source, target):
    """
    """

    path = [target]
    #print 'initial path: ', path

    while path[-1] != source:
        #print 'iterated path: ', path[-1]
        path.append(prev[path[-1]])
        #print 'append: ', prev[path[-1]]

    path.reverse()
    return path
