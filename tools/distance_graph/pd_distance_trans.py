#!/usr/bin/env python

from libpl.pdstor import PDStoreExpander

from libpl.pdstor import ClassifyDatabase

import timeit

import cPickle as pickle
import networkx as nx
import os

DATADIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                       "plCurve", "data", "pdstors")
#DATADIR = "/home/harrison/Projects/plCurve/data/pdstors"

def homfly_graph(num_crossings):
    G = nx.Graph()

    for pd_graph in shadow_graph_iter(num_crossings):
        for uid_a, uid_b, data in pd_graph.edges(data=True):
            pd_a, pd_b = data['pds']
            hf_a, hf_b = pd_a.homfly(), pd_b.homfly()

            G.add_edge(hf_a, hf_b,
                       weight=G.get_edge_data(hf_a,hf_b,{'weight':0})['weight']+data['weight'])
        del pd_graph

    return G

def shadow_graph_iter(num_crossings):
    thin = True
    expander = PDStoreExpander(dirloc=DATADIR, debug=False)
    pdstor = expander.open_shadows(num_crossings, thin=thin, debug=True,
                                   orient_all=False, homflys=False, max_components=1)
    for i, shadow in enumerate(pdstor):
        G = nx.Graph()
        pds = dict()
        for pdcode, homfly, uid in shadow:
            if pdcode.ncomps > 1:
                continue # We only want knots for now
            G.add_node(uid,pd=pdcode)
            pds[uid[2]] = pdcode
        print "Diagram distance graph %s loaded with %s nodes"%(i,G.number_of_nodes())
        for start_uid, data in G.nodes(data=True):
            start_pd = data['pd']
            ncross, pos, crs_sign, cmp_sign = start_uid
            for idx in range(ncross):
                new_x = tuple((x+1)%2 if i==idx else x for i,x in enumerate(crs_sign))
                end_uid = (ncross, pos, new_x, cmp_sign)
                end_pd = pds[new_x]
                G.add_edge(start_uid, end_uid, pds=(start_pd, end_pd),
                           weight=G.get_edge_data(start_uid, end_uid,{'weight':0})['weight']+1)
        print "Edges of diagram distance %s graph computed: %s edges"%(i,G.number_of_edges())
        yield G

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
