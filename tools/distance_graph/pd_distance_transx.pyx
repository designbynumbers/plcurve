# cython: profile=False
# filename: pd_distance_transx.pyx

from libpl.pdstor import PDStoreExpander
from libpl.pdstor import ClassifyDatabase

from libpl.pdcode import HOMFLYPolynomial

import timeit
from cpython cimport bool

import cPickle as pickle
import networkx as nx
import os

DATADIR = os.path.join("plCurve", "data", "pdstors")
#DATADIR = "/home/harrison/Projects/plCurve/data/pdstors"

COLLISION_HOMFLY = (
    HOMFLYPolynomial("-a^{-8} + -2a^{-6} + 2a^{-6}z^{2} + -a^{-4} + a^{-4}z^{2} + -a^{-4}z^{4} + -a^{-2}z^{2} + a^{-2}z^{4} + 1 + -z^{2}"),
    HOMFLYPolynomial("1 + -z^{2} + -a^{2}z^{2} + a^{2}z^{4} + -a^{4} + a^{4}z^{2} + -a^{4}z^{4} + -2a^{6} + 2a^{6}z^{2} + -a^{8}"),
)

cpdef object homfly_graph(object num_crossings):
    cdef object G = nx.Graph()
    cdef bool thin = False
    cdef object expander = PDStoreExpander(dirloc=DATADIR, debug=False)
    cdef object pdstor = expander.open_shadows(num_crossings, thin=thin, debug=True,
                                   orient_all=False, homflys=False, max_components=1)
    #i = 0
    cdef long i = 0
    cdef long j,k, idx
    cdef long ncross, pos

    cdef object pd_a, pd_b
    cdef tuple uid, end_uid, start_uid, crs_sign, cmp_sign
    cdef dict pds, data

    #cdef object shadow
    #for i, shadow in enumerate(pdstor):
    for shadow in pdstor:
        pd_graph = nx.Graph()
        pds = dict()
        for pdcode, homfly, uid in shadow:
            if pdcode.ncomps > 1:
                continue # We only want knots for now
            pd_graph.add_node(uid,pd=pdcode)
            pds[uid[2]] = pdcode
        #print "Diagram distance graph %s loaded with %s nodes"%(i,pd_graph.number_of_nodes())
        for start_uid, data in pd_graph.nodes(data=True):
            start_pd = data['pd']
            ncross, pos, crs_sign, cmp_sign = start_uid
            for idx in range(ncross):
                new_x = list() #tuple((x+1)%2 if i==idx else x for i,x in enumerate(crs_sign))
                for k,x in enumerate(crs_sign):
                    new_x.append( (x+1)%2 if k==idx else x )
                new_x = tuple(new_x)

                end_uid = (ncross, pos, new_x, cmp_sign)
                end_pd = pds[new_x]
                try:
                    pd_graph[start_uid][end_uid]['weight'] += 1
                except KeyError:
                    pd_graph.add_edge(start_uid, end_uid, pds=(start_pd, end_pd), weight=1)
        #print "Edges of diagram distance %s graph computed: %s edges"%(i,pd_graph.number_of_edges())

        for uid_a, uid_b, data in pd_graph.edges(data=True):
            pd_a, pd_b = data['pds']
            hf_a, hf_b = pd_a.homfly(), pd_b.homfly()

            if hf_a in COLLISION_HOMFLY:
                key_a = (hf_a, len(pd_a.snappy_manifold().splitting_surfaces()))
            else:
                key_a = hf_a

            if hf_b in COLLISION_HOMFLY:
                key_b = (hf_b, len(pd_b.snappy_manifold().splitting_surfaces()))
            else:
                key_b = hf_b

            try:
                G[key_a][key_b]['weight'] += data['weight']
            except KeyError:
                G.add_edge(key_a, key_b, weight=data['weight'])

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
    for specA, specB, edata in hG.edges(data=True):
        if isinstance(specA, HOMFLYPolynomial):
            homA = specA
            resA = cdb.classify_homfly(homA, N)
            nodeA = resA[0].concise()
        else:
            if specA[0] == COLLISION_HOMFLY[0]:
                if specA[1] == 2:
                    nodeA = "4_1*#5_2"
                else:
                    nodeA = "9_12"
            else:
                if specA[1] == 2:
                    nodeA = "4_1#5_2*"
                else:
                    nodeA = "9_12*"

        if isinstance(specB, HOMFLYPolynomial):
            homB = specB
            resB = cdb.classify_homfly(homB, N)
            nodeB = resB[0].concise()
        else:
            if specB[0] == COLLISION_HOMFLY[0]:
                if specB[1] == 2:
                    nodeB = "4_1*#5_2"
                else:
                    nodeB = "9_12"
            else:
                if specB[1] == 2:
                    nodeB = "4_1#5_2*"
                else:
                    nodeB = "9_12*"

        G.add_edge(nodeA, nodeB, weight=edata['weight'])
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
