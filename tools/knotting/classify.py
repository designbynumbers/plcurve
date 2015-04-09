from database import *
from itertools import izip
from collections import defaultdict
import gzip

def unmattify(mktype):
    kt = mktype.strip().replace("m","*")
    if "{" in kt:
        return kt.replace("{","").replace(", ","#").replace("\"","").replace("}","").replace(",",", ")
    else:
        return kt

if __name__ == "__main__":
    session = Session()
    N = 7
    with open("../../data/KnotTheorycodes/%s.dat"%N) as f:
        matt_types = [unmattify(line) for line in f]
    with gzip.open("../../data/KnotTheorycodes/%s.knottheory.gz"%N) as f:
        matt_2pds = [line.strip() for line in f]
    matt_pds = [matt_2pds[2*i:2*i+2] for i in range(len(matt_2pds)/2)]
    #print matt_pds[0]

    db_diagrams = session.query(Diagram).join(Shadow).filter(Shadow.n_cross==N).all()
    knot_uids = sorted([uid for uid, in session.query(Shadow.uid).filter(Shadow.n_cross==N).all()])
    #print knot_uids

    nknotshadow = len(matt_types)/(2**N)/2
    #print nknotshadow
    #print "----"
    matt_dict = [[] for _ in range(nknotshadow)]
    for i,uid_list in enumerate(matt_dict):
        for j in range(2):
            xmasked = []
            for k in range(2**N):
                iota = k+2**N*j+2**(N+1)*i
                xmasked.append((matt_types[iota],)+tuple(matt_pds[iota]))
            uid_list.append(xmasked)

    #print matt_dict
#    sorted_diagrams = sorted(db_diagrams, key=lambda d: 
    
    for db_diagram in db_diagrams:
        #print db_diagram.shadow.uid
        #print db_diagram.cross_mask, 2**N-1-db_diagram.cross_mask
        #print db_diagram.comp_mask
        factorizations = session.query(LinkFactorization)\
                     .filter(LinkFactorization.n_cross<=db_diagram.shadow.n_cross)\
                     .filter(LinkFactorization.homfly==db_diagram.homfly).all()
        assert len(factorizations) == 1
        fzn, = factorizations

        mtype, comment, pdcode = matt_dict[knot_uids.index(db_diagram.shadow.uid)][db_diagram.comp_mask][2**N-1-db_diagram.cross_mask]
        #matt_types[knot_uids.index(db_diagram.shadow.uid)*(2**N+1)+
        #           db_diagram.comp_mask*(2**N)+(2**N-db_diagram.cross_mask-1)]
                             
        if str(fzn) != mtype:
            print fzn, mtype
            print comment
            print pdcode
            print db_diagram.shadow.uid, db_diagram.cross_mask_list
            #print matt_dict[db_diagram.shadow.uid]
            print
        pass
