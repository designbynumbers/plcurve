from database import *
from libpl.pdstor import *
from collections import defaultdict
P_unk = HOMFLYPolynomial('1')

def load_diagrams(session, n_cross):
    print "There are %s shadows to process"%session.query(Shadow)\
                                  .filter(Shadow.n_cross==n_cross).count()
    for i, db_shadow in enumerate(session.query(Shadow)\
                                  .filter(Shadow.n_cross==n_cross)):
        if i%100 == 0:
            print i
            session.commit()

        # Two knot diagrams can only be isotopic if they have the same (absolute) total sign.
        iso_classes = defaultdict(dict)
        for comp_pd, cmask in PDStoreExpander.component_combinations(db_shadow.pd, True, thin=False):
            for pd, xmask in PDStoreExpander.crossing_combinations(comp_pd):
                sign_count = sum(x*2-1 for x in xmask)
                if not (sign_count in iso_classes and pd in iso_classes[sign_count]):
                    homfly = pd.homfly()

                    db_iso = DiagramClass(
                        shadow = db_shadow,
                        homfly = homfly,
                    )
                    session.add(db_iso)
                    iso_classes[sign_count][pd] = db_iso
                else:
                    db_iso = iso_classes[sign_count][pd]

                db_pd = Diagram(
                    iso_class = db_iso,
                    cross_mask = xmask,
                    comp_mask = cmask,
                    )
                session.add(db_pd)


    session.commit()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Expand shadows in the database into diagrams.")
    parser.add_argument('crossings', metavar='N', type=int, nargs='+',
                        help="crossing counts shadows to expand")
    args = parser.parse_args()

    session = Session()
    for n_cross in args.crossings:
        print n_cross
        load_diagrams(session, n_cross)
    session.close()
