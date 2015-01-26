# cython: profile=True
# filename: pd_distance_transx.pyx

from libpl.pdstor import PDStoreExpander

from libpl.pdstor import ClassifyDatabase

import timeit
from cpython cimport bool

import cPickle as pickle
import networkx as nx
import igraph as ig
import os

DATADIR = os.path.join("plCurve", "data", "pdstors")
#DATADIR = "/home/harrison/Projects/plCurve/data/pdstors"

cpdef object homfly_graph(object num_crossings):
    cdef object G = ig.Graph()
    cdef bool thin = True
    cdef object expander = PDStoreExpander(dirloc=DATADIR, debug=False)
    cdef object pdstor = expander.open_shadows(num_crossings, thin=thin, debug=True,
                                   orient_all=False, homflys=False, max_components=1)
    #cdef int i
    #cdef object shadow
    for i, shadow in enumerate(pdstor):
        pd_graph = ig.Graph()
        pds = dict()

        uid_verts = dict()
        for pdcode, homfly, uid in shadow:
            if pdcode.ncomps > 1:
                continue # We only want knots for now
            uid_verts[uid] = pd_graph.vcount()
            pd_graph.add_vertex(name=uid,pd=pdcode)
            pds[uid[2]] = pdcode
        print "Diagram distance graph %s loaded with %s nodes"%(i,pd_graph.vcount())

        for start_v in pd_graph.vs:
            start_pd = start_v['pd']
            start_uid = start_v['name']
            ncross, pos, crs_sign, cmp_sign = start_uid
            for idx in range(ncross):
                new_x = list() #tuple((x+1)%2 if i==idx else x for i,x in enumerate(crs_sign))
                for j,x in enumerate(crs_sign):
                    new_x.append( (x+1)%2 if j==idx else x )
                new_x = tuple(new_x)
                end_uid = (ncross, pos, new_x, cmp_sign)
                end_pd = pds[new_x]

                start_vid, end_vid = start_v.index, uid_verts[end_uid]
                k = pd_graph.get_eid(start_vid, end_vid, error=False)
                if k == -1:
                    pd_graph.add_edge(start_vid, end_vid,
                                      pds=(start_pd, end_pd),
                                      weight=1)
                else:
                    pd_graph.es[k]['weight'] += 1

        print "Edges of diagram distance %s graph computed: %s edges"%(i,pd_graph.ecount())

        verts = dict()
        for edge in pd_graph.es:
            pd_a, pd_b = edge['pds']
            hf_a, hf_b = pd_a.homfly(), pd_b.homfly()

            if hf_a not in verts:
                vid_a = G.vcount()
                v_a = G.add_vertex(name=hf_a)
                verts[hf_a] = vid_a
            else:
                vid_a = verts[hf_a]

            if hf_b not in verts:
                vid_b = G.vcount()
                v_b = G.add_vertex(name=hf_b)
                verts[hf_b] = vid_b
            else:
                vid_b = verts[hf_b]


            k = G.get_eid(vid_a, vid_b, error=False)
            if k == -1:
                G.add_edge(vid_a, vid_b, weight=edge['weight'])
            else:
                G.es[k]['weight'] += edge['weight']

        del pd_graph

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
