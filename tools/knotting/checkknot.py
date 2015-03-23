from database import *
from libpl.pdstor import *
P_unk = HOMFLYPolynomial('1')

if __name__ == "__main__":
    session = Session()
    n_cross = 6
    for db_shadow in session.query(Shadow).filter(Shadow.n_cross==n_cross):
        for pd, xmask in PDStoreExpander.crossing_combinations(db_shadow.pd):
            homfly = pd.homfly()
            if homfly != P_unk:
                (ret, ), = session.query(exists()
                .where(Diagram.shadow==db_shadow)
                .where(Diagram.cross_mask==xmask))
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
                    print pd.crossings, pd.homfly()
    session.commit()
    session.close()
