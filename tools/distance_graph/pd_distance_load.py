import networkx as nx
from libpl.pdcode import *
import cPickle as pickle
from libpl.pdstor import ClassifyDatabase

from pd_distance import draw_graph_subgraph, draw_graph, homfly_graph_to_knots


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Count isotopy classes of diagrams from shadows")
    parser.add_argument('num_crossings', nargs='+', metavar='N', type=int,
                        help="the crossing number of the shadows file to use")
    args = parser.parse_args()
    #graph = populate_graph(args.num_crossings, session)
    #graph2 = populate_graph_uid(args.num_crossings, session)
    #draw_graph(graph)

    n, k = args.num_crossings


    print "Loading classification database.."
    cdb = ClassifyDatabase()
    cdb.load(load_links=False)
    cdb.calculate_composites(n, num_components=1)
    print "Done! Computing graph..."

    with open("%sx_hfly_dgraph.pickle"%(n,), "r") as f:
        Hn = pickle.Unpickler(f).load()
    with open("%sx_hfly_dgraph.pickle"%(k,), "r") as f:
        Hk = pickle.Unpickler(f).load()

    Gn = homfly_graph_to_knots(Hn, cdb, n)
    Gk = homfly_graph_to_knots(Hk, cdb, k)

    print "Graph loaded: %s edges"%len(Gn.edges())

    ambiguous_nodes = ["Knot[9, 12]", "Knot[9, 12]*"]
    Gn.remove_nodes_from(ambiguous_nodes)

    full_paths = nx.shortest_path_length(Gn)
    sub_paths = nx.shortest_path_length(Gk)

    subgraph_paths = dict()
    for source in Gk:
        subgraph_paths[source] = dict()
        for target in Gk:
            subgraph_paths[source][target] = nx.shortest_path_length(
                Gn, source=source, target=target)

    for source, targets in subgraph_paths.iteritems():
        for target, length in targets.iteritems():
            longer = sub_paths[source][target]
            #print source, target, length, longer
            assert length <= longer
            if length < longer:
                print "{:>30} --> {:<30}.. Was {:2}, now {:2}".format(
                    source, target, longer, length)
                print nx.shortest_path(Gn, source, target)
                print

    for homfly in Hn:
        res = cdb.classify_homfly(homfly, n)
        if len(res) != 1:
            print homfly
            for kt in res:
                print kt
            print "====="

    draw_graph_subgraph(Gn, Gk)
    #print full_paths["Knot[6, 1]"]

    # TODO:
    # Distances to unknot
    # KnotInfo
    # Find out if 8x and 9x graphs show known unknotting distances for (8x, 9x) knots
    # KnotPlot... has knot distance data built in
    # http://homepage.math.uiowa.edu/~idarcy/TAB/tabnov.pdf
    # Try to compare graph results to this data (check some unknown values)
    # Set a machine to 10.pdstor?
