from database import *
from collections import defaultdict

def ktname_to_header(ktname):
    # Turns Knot[4, 1] into 4.1 and Knot[5, 1]* into 5.1m
    return ktname.replace("Knot[","")\
                 .replace("]","")\
                 .replace(", ",".")\
                 .replace("*","m")\
                 .replace("#","c")

def get_monobi_dict(session, n_cross):
    total = session.query(Shadow).filter(Shadow.n_cross==n_cross).count()
    monos = [shad.pd.faces[shad.pd.nfaces-1].nedges==1 for
             shad in session.query(Shadow).filter(Shadow.n_cross==n_cross)].count(True)
    bigos = [any(f.nedges==2 for f in shad.pd.faces) for
             shad in session.query(Shadow).filter(Shadow.n_cross==n_cross)].count(True)
    mobis = [(any(f.nedges==2 for f in shad.pd.faces) and
              shad.pd.faces[shad.pd.nfaces-1].nedges==1) for
             shad in session.query(Shadow).filter(Shadow.n_cross==n_cross)].count(True)
    moobis = [(any(f.nedges==2 for f in shad.pd.faces) or
              shad.pd.faces[shad.pd.nfaces-1].nedges==1) for
             shad in session.query(Shadow).filter(Shadow.n_cross==n_cross)].count(True)
    print total, monos, bigos, mobis, moobis

from sqlalchemy import func

def get_mono_dict(session, n_cross):
    n_shadows = session.query(Shadow).filter(Shadow.n_cross==n_cross).count()
    n_diagrams = session.query(DiagramClass).join(Shadow).filter(Shadow.n_cross==n_cross).count()
    diagrams = defaultdict(lambda: 0.0)
    shadows = defaultdict(lambda: 0.0)
    #print "about to check"
    for pd, count in session.query(Shadow.pd, func.count(DiagramClass.id)).join(Shadow.diagram_classes).filter(Shadow.n_cross==n_cross).group_by(Shadow.id):
        N = [len(face) for face in pd.faces].count(1)
        shadows[N] += 1
        diagrams[N] += count

    print str(n_cross) + "\t" + "\t".join(str(diagrams[i]/n_diagrams) for i in range(0,10))

def get_bigo_dict(session, n_cross):
    n_shadows = session.query(Shadow).filter(Shadow.n_cross==n_cross).count()
    n_diagrams = session.query(DiagramClass).join(Shadow).filter(Shadow.n_cross==n_cross).count()
    diagrams = defaultdict(lambda: 0.0)
    shadows = defaultdict(lambda: 0.0)
    for pd, count in session.query(Shadow.pd, func.count(DiagramClass.id)).join(Shadow.diagram_classes).filter(Shadow.n_cross==n_cross).group_by(Shadow.id):
        N = [len(face) for face in pd.faces].count(2)
        shadows[N] += 1
        diagrams[N] += count

    print str(n_cross) + "\t" + "\t".join(str(diagrams[i]/n_diagrams) for i in range(0,10))

def reduce_pd(pdcode):
    loops = list(pdcode.get_R1_loops())
    if not loops:
        return pdcode
    return reduce_pd(pdcode.R1_loop_deletion(loops[0].vertices[0]))

def count_treelike(session, n_cross):
    n_shadows = session.query(Shadow).filter(Shadow.n_cross==n_cross).count()
    n_diagrams = session.query(DiagramClass).join(Shadow).filter(Shadow.n_cross==n_cross).count()
    diagrams = 0.0
    shadows = 0.0
    for pd, count in session.query(Shadow.pd, func.count(DiagramClass.id)).join(Shadow.diagram_classes).filter(Shadow.n_cross==n_cross).group_by(Shadow.id):
        pd.set_all_crossing_signs((0,)*n_cross)
        if reduce_pd(pd).ncross == 0:
            shadows += 1
            diagrams += count

    print str(n_cross) + "\t" + str(shadows/n_shadows) + "\t" + (str(diagrams/n_diagrams))

def avg_reduced_cn(session, n_cross):
    n_shadows = session.query(Shadow).filter(Shadow.n_cross==n_cross).count()
    n_diagrams = session.query(DiagramClass).join(Shadow).filter(Shadow.n_cross==n_cross).count()
    diagrams = defaultdict(lambda: 0.0)
    shadows = defaultdict(lambda: 0.0)
    for pd, count in session.query(Shadow.pd, func.count(DiagramClass.id)).join(Shadow.diagram_classes).filter(Shadow.n_cross==n_cross).group_by(Shadow.id):
        pd.set_all_crossing_signs((0,)*n_cross)
        N = reduce_pd(pd).ncross
        shadows[N] += 1
        diagrams[N] += count

    shad_avg = sum(N*count for N, count in shadows.iteritems())/n_shadows
    dias_avg = sum(N*count for N, count in diagrams.iteritems())/n_diagrams
    print str(n_cross) + "\t" + str(shad_avg) + "\t" + "\t".join(str(shadows[i]/n_shadows) for i in range(0,10))
    print str(n_cross) + "\t" + str(dias_avg) + "\t" + "\t".join(str(diagrams[i]/n_diagrams) for i in range(0,10))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Count knots with monogons/bigons")
    parser.add_argument('max_crossings', type=int,
                        help="Largest number of crossings to analyze")
    args = parser.parse_args()

    session = Session()
    mbdicts = {}
    """
    print "n\ttreeshad\ttreedia"
    for n_cross in range(3, args.max_crossings+1):
        count_treelike(session, n_cross)

    print "n\t" + "\t".join(str(k) for k in range(0,10))
    for n_cross in range(3, args.max_crossings+1):
        mbdicts[n_cross] = get_mono_dict(session, n_cross)

    print "n\t" + "\t".join(str(k) for k in range(0,10))
    for n_cross in range(3, args.max_crossings+1):
        mbdicts[n_cross] = get_bigo_dict(session, n_cross)
    """
    print "n\tavg\t" + "\t".join(str(k) for k in range(0,10))
    for n_cross in range(3, args.max_crossings+1):
        avg_reduced_cn(session, n_cross)
