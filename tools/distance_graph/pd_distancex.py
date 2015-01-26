CYTHON = True
if CYTHON:
    import pyximport
    pyximport.install()
    from pd_distance_transx import *
else:
    from pd_distance_trans import *

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Count isotopy classes of diagrams from shadows")
    parser.add_argument('num_crossings', nargs='+', metavar='N', type=int,
                        help="the crossing number of the shadows file to use")
    args = parser.parse_args()

    n = args.num_crossings

    #print "Loading classification database.."
    #cdb = ClassifyDatabase()
    #cdb.load()
    #cdb.calculate_composites(max(n))
    #print "Done! Computing graph..."

    hfly_graph = homfly_graph(n)
    #knot_graph = homfly_graph_to_knots(hfly_graph, cdb, max(n))

    #print [(u,v,edata['weight']) for u,v,edata in knot_graph.edges(data=True)]

    with open("%sx_hfly_dgraph_transitions.pickle"%("_".join(str(x) for x in n)), "wb") as f:
        f.write(pickle.dumps(hfly_graph, protocol=pickle.HIGHEST_PROTOCOL))
