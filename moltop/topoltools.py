import numpy as np
import matplotlib.pyplot as plt

def harmonic_energy(pos, graph):

    pos = pos.reshape((len(graph), 3))

    en = 0.
    for vertex in graph:
        for edge, dist in zip(graph[vertex], graph.dists[vertex]):
            en += (np.linalg.norm(pos[vertex] - pos[edge]) - dist)**2

        #for vx in graph:
        #    #this looks like a funky pentagram
        #    #en += np.linalg.norm(pos[vertex] - pos[vx])**2
        #    #
        #    if vx not in graph.edges[vertex]:
        #        #the distance is not a good measure here, get actual distance by calculating path length
        #        en += (np.linalg.norm(pos[vertex] - pos[edge]) - 1.)**2

    return en

def plot_graph(pos, graph):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    print pos[0]
    for vertex, crd in zip(graph, pos):
        ax.scatter(crd[0], crd[1], crd[2], s=20, marker='$%i$' %vertex, color='r')

    #plt.plot(pos[:,0], pos[:,1], 'o')

    for vertex in graph:
        for edge in graph[vertex]:
            ln = np.vstack([pos[vertex], pos[edge]])
            ax.plot(ln[:,0], ln[:,1], ln[:,2], color='k')

