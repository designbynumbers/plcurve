from database import *
from libpl.pdstor import *
P_unk = HOMFLYPolynomial('1')

def load_diagrams(session, n_cross):
    for i, db_shadow in enumerate(session.query(Shadow).filter(Shadow.n_cross==n_cross)):
        if i%100 == 0:
            print i
            session.commit()

        iso_classes = []

        # Two knot diagrams can only be isotopic if they have the same (absolute) total sign.
        sign_counts = []
        for pd, xmask in PDStoreExpander.crossing_combinations(db_shadow.pd):
            sign_count = abs(sum(x*2-1 for x in xmask))
            if (sign_count in sign_counts and
                any((pd.isotopic(iso_class) for iso_class in iso_classes))):
                ## Still have to add diagram :(
                continue
            else:
                sign_counts.append(sign_count)
                iso_classes.append(pd)
            
            homfly = pd.homfly()
            ret = False
            #(ret, ), = session.query(exists()
            #.where(Diagram.shadow==db_shadow)
            #.where(Diagram.cross_mask==xmask))
            #ret = session.query(Diagram)\
            #             .filter(Diagram.shadow==db_shadow)\
            #             .filter(Diagram.cross_mask==xmask).all()
            #print ret
            if not ret:
                iso_pd = DiagramClass(
                    shadow = db_shadow,
                    homfly = homfly,
                )
                session.add(iso_pd)
                db_pd = Diagram(
                    iso_class = iso_pd,
                    cross_mask = xmask,
                    comp_mask = (0,)
                )
                session.add(db_pd)
                #print pd.crossings, pd.homfly()
    session.commit()

if __name__ == "__main__":
    session = Session()
    for n_cross in range(3,6+1):
        print n_cross
        load_diagrams(session, n_cross)
    session.close()
