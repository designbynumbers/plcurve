import pstats, cProfile

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

    cProfile.runctx("homfly_graph(n)", globals(), locals(), "Profile.prof")

    s = pstats.Stats("Profile.prof")
    s.strip_dirs().sort_stats("time").print_stats()
