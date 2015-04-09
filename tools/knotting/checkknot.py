from database import *
from libpl.pdstor import *
P_unk = HOMFLYPolynomial('1')

def load_diagrams(session, n_cross):
    for i, db_shadow in enumerate(session.query(Shadow).filter(Shadow.n_cross==n_cross)):
        if i%100 == 0:
            print i
        for pd, xmask in PDStoreExpander.crossing_combinations(db_shadow.pd):
            homfly = pd.homfly()
            if homfly != P_unk:
                ret = False
                #(ret, ), = session.query(exists()
                #.where(Diagram.shadow==db_shadow)
                #.where(Diagram.cross_mask==xmask))
                #ret = session.query(Diagram)\
                #             .filter(Diagram.shadow==db_shadow)\
                #             .filter(Diagram.cross_mask==xmask).all()
                #print ret
                if not ret:
                    db_pd = Diagram(
                        shadow = db_shadow,
                        homfly = homfly,
                        cross_mask = xmask,
                        comp_mask = (0,)
                    )
                    session.add(db_pd)
                    #print pd.crossings, pd.homfly()
    session.commit()

if __name__ == "__main__":
    session = Session()
    for n_cross in range(8,9):
        print n_cross
        load_diagrams(session, n_cross)
    session.close()
