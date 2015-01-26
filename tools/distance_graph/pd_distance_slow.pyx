# cython: profile=True
# filename: pd_distance_transx.pyx

from libpl.pdstor import PDStoreExpander

from libpl.pdstor import ClassifyDatabase

import timeit
from cpython cimport bool

import cPickle as pickle
import networkx as nx
import os

DATADIR = os.path.join("plCurve", "data", "pdstors")
#DATADIR = "/home/harrison/Projects/plCurve/data/pdstors"

cpdef object homfly_graph(object num_crossings):
    cdef object G = nx.Graph()
    cdef bool thin = True
    cdef object expander = PDStoreExpander(dirloc=DATADIR, debug=False)
    cdef object pdstor = expander.open_shadows(num_crossings, thin=thin, debug=True,
                                   orient_all=False, homflys=False, max_components=1)
    #i = 0
    cdef long i = 0
    cdef long j,k
    cdef long ncross, pos
    #cdef object shadow
    #for i, shadow in enumerate(pdstor):
    for shadow in pdstor:
        pd_graph = nx.Graph()
        pds = dict()

        for start_pd, homfly, start_uid in shadow:
            pds[start_uid[2]] = start_pd

            ncross, pos, crs_sign, cmp_sign = start_uid
            for idx in range(ncross):
                new_x = list()
                for k,x in enumerate(crs_sign):
                    new_x.append( (x+1)%2 if k==idx else x )
                new_x = tuple(new_x)

                end_pd = start_pd.copy()
                end_pd.toggle_crossing(idx)
                #end_pd = PDStoreExpander.crossing_sign_mask(start_pd, crs_sign)
                end_uid = (ncross, pos, new_x, cmp_sign)

                try:
                    pd_graph[start_uid][end_uid]['weight'] += 1
                except KeyError:
                    pd_graph.add_edge(start_uid, end_uid, pds=(start_pd, end_pd), weight=1)

        for uid_a, uid_b, data in pd_graph.edges(data=True):
            pd_a, pd_b = data['pds']
            hf_a, hf_b = pd_a.homfly(), pd_b.homfly()

            try:
                G[hf_a][hf_b]['weight'] += data['weight']
            except KeyError:
                G.add_edge(hf_a, hf_b, weight=data['weight'])

        del pd_graph
        i += 1

    return G

def draw_graph(G):
    import matplotlib.pyplot as plt
    nx.draw_networkx(G, with_labels=True)

    plt.show()

    #import IPython
    #IPython.embed()

def draw_graph_subgraph(G, H):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    nodes = G.nodes()
    node_prefs = list((800, "b") if x in H else (300, "r") for x in nodes)
    node_sizes, node_colors = zip(*node_prefs)

    edges = G.edges()
    edge_styles = list("b" if H.has_edge(*e) else "r" for
                       e in edges)

    nx.draw_networkx(G, with_labels=True, nodelist=nodes, edgeslist=edges,
                     node_size=node_sizes,
                     node_color=node_colors,
                     edge_color=edge_styles,
    )

    plt.subplots_adjust(top=1,left=0,right=1,bottom=0)
    #fig.tight_layout()
    #ax.axis('off')
    plt.show()

def homfly_graph_to_knots(hG, cdb, N):
    G = nx.Graph()
    for homA, homB, edata in hG.edges(data=True):
        resA = cdb.classify_homfly(homA, N)
        resB = cdb.classify_homfly(homB, N)
        G.add_edge(resA[0].concise(), resB[0].concise(), weight=edata['weight'])
    return G

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Count isotopy classes of diagrams from shadows")
    parser.add_argument('num_crossings', nargs='+', metavar='N', type=int,
                        help="the crossing number of the shadows file to use")
    args = parser.parse_args()

    n = args.num_crossings

    print "Loading classification database.."
    cdb = ClassifyDatabase()
    cdb.load()
    cdb.calculate_composites(max(n))
    print "Done! Computing graph..."

    hfly_graph = homfly_graph(n)
    knot_graph = homfly_graph_to_knots(hfly_graph, cdb, max(n))

    #print [(u,v,edata['weight']) for u,v,edata in knot_graph.edges(data=True)]

    with open("%sx_hfly_dgraph_transitions.pickle"%("_".join(str(x) for x in n)), "wb") as f:
        f.write(pickle.dumps(hfly_graph, protocol=pickle.HIGHEST_PROTOCOL))
