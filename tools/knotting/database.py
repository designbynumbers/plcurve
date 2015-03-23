from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref, column_property
import sqlalchemy.types as types
from sqlalchemy.ext.hybrid import hybrid_property
from itertools import compress

engine = create_engine("sqlite:///shadows.db", echo=True)

from libpl.pdcode import *
class SerializedDiagram(types.TypeDecorator):
    impl = types.LargeBinary

    def process_bind_param(self, value, dialect):
        if isinstance(value, PlanarDiagram):
            return value.serialize()
        return value

    def process_result_value(self, value, dialect):
        if value is not None:
            return PlanarDiagram.deserialize(value)
        return value

def bin_list_to_int(blist):
    N = len(blist)-1
    return sum(n*(2**(N-i)) for i,n in enumerate(blist))
def int_to_bin_list(bint, n=None):
    ret = tuple(int(x) for x in bin(bint)[2:])
    if n is None:
        return ret
    return (0,)*(n-len(ret)) + ret

class Mask(types.TypeDecorator):
    impl = types.Integer

    def process_bind_param(self, value, dialect):
        if isinstance(value, list) or isinstance(value, tuple):
            return bin_list_to_int(value)
        return value

class HOMFLYType(types.TypeDecorator):
    impl = types.String

    def process_bind_param(self, value, dialect):
        if value is not None:
            value = str(value)
        return value

    def process_result_value(self, value, dialect):
        if value is not None:
            value = HOMFLYPolynomial(value)
        return value

Base = declarative_base()
class Shadow(Base):
    __tablename__ = 'shadows'

    id = Column(Integer, primary_key=True)
    uid = Column(Integer)

    n_cross = Column(Integer)
    n_comps = Column(Integer)

    hash = Column(String)

    pd = Column(SerializedDiagram)

    def __repr__(self):
        return "<Shadow(n_cross=%s, n_comps=%s, uid=%s, hash='%s')>"%(
            self.n_cross, self.n_comps, self.uid, self.hash)

class Diagram(Base):
    __tablename__ = 'diagrams'

    id = Column(Integer, primary_key=True)
    shadow_id = Column(Integer, ForeignKey('shadows.id'))
    shadow = relationship("Shadow", backref=backref('diagrams', order_by=id))

    n_cross = column_property(shadow.n_cross)

    homfly = Column(HOMFLYType)

    cross_mask = Column(Mask)
    comp_mask = Column(Mask)

    @property
    def cross_mask_list(self):
        return int_to_bin_list(self.cross_mask, self.shadow.n_cross)
    @property
    def comp_mask_list(self):
        return int_to_bin_list(self.comp_mask, self.shadow.n_comps)

    @property
    def pd(self):
        pd = self.shadow.pd
        unsignature = (1-x for x in self.comp_mask_list)
        for comp in compress(range(pd.ncomps), unsignature):
           pd.reorient_component(comp, 0)
        pd.set_all_crossing_signs(self.cross_mask_list)
        return pd

    def __repr__(self):
        return "<Diagram(n_cross=%s, n_comps=%s, shadow=%s, cross_mask=%s, comp_mask=%s)>"%(
            self.shadow.n_cross, self.shadow.n_comps, self.shadow.uid,
            self.cross_mask, self.comp_mask)

class KnotFactor(Base):
    __tablename__ = 'factors'

    id = Column(Integer, primary_key=True)
    homfly = Column(HOMFLYType)

    n_cross = Column(Integer)
    n_comps = Column(Integer)

Base.metadata.create_all(engine)
Session = sessionmaker(bind=engine)

if __name__ == "__main__":
    session = Session()

    from libpl.pdstor import *

    """
    expander = PDStoreExpander(dirloc="../../data/pdstors", debug=False)
    pdstor = expander.open_shadows([4], debug=True, orient_all=False, homflys=False, max_components=1)
    for i,shadow in enumerate(pdstor):
        db_shadow = None
        for pdcode, homfly, uid in shadow:
            if db_shadow is None:
                db_shadow = Shadow(
                    uid=pdcode.uid,
                    hash=pdcode.hash,
                    n_cross=pdcode.ncross,
                    n_comps=pdcode.ncomps,
                    pd=pdcode)
                session.add(db_shadow)
            db_pdcode = Diagram(
                homfly=str(pdcode.homfly()),
                shadow=db_shadow,
                cross_mask=uid[2],
                comp_mask=uid[3])
            session.add(db_pdcode)
    """

    N = 9
    with open("../../data/pdstors/%s.pdstor"%N, "rb") as f:
        for shadow in PlanarDiagram.read_all(f, read_header=True):
            if shadow.ncomps > 1:
                continue
            (ret, ), = session.query(exists().where(
                Shadow.uid==shadow.uid).where(
                    Shadow.n_cross==shadow.ncross))
            if ret:
                continue
            db_shadow = Shadow(
                uid = shadow.uid,
                hash = shadow.hash,
                n_cross = shadow.ncross,
                n_comps = shadow.ncomps,
                pd = shadow
            )
            session.add(db_shadow)

    session.commit()
    session.close()
