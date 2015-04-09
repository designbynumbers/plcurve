from database import *

if __name__ == "__main__":
    session = Session()
    for db_diagram in session.query(Diagram).join(Shadow).filter(Shadow.n_cross==6):
        factorizations = session.query(LinkFactorization)\
                     .filter(LinkFactorization.n_cross<=db_diagram.shadow.n_cross)\
                     .filter(LinkFactorization.homfly==db_diagram.homfly).all()
        
        print [str(fzn) for fzn in factorizations]
        pass
