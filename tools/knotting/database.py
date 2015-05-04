from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref, column_property
import sqlalchemy.types as types
from sqlalchemy.ext.hybrid import hybrid_property
from itertools import compress

engine = create_engine("sqlite:///shadows.db", echo=False)

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

class FactorizationFactor(Base):
    __tablename__ = 'factorization_factors'
    product_id = Column(Integer, ForeignKey('factorizations.id'), primary_key=True)
    factor_id = Column(Integer, ForeignKey('factors.id'), primary_key=True)
    factor = relationship("LinkFactor")
    multiplicity = Column(Integer)

diagram_factorizations = Table(
    'diagram_factorizations', Base.metadata,
    Column('diagram_id', Integer, ForeignKey('diagrams.id')),
    Column('factorization_id', Integer, ForeignKey('factorizations.id'))
)

class Shadow(Base):
    __tablename__ = 'shadows'

    id = Column(Integer, primary_key=True)
    uid = Column(Integer)

    n_cross = Column(Integer, index=True) # # crossings
    n_comps = Column(Integer, index=True) # # components
    n_autos = Column(Integer) # # automorphisms

    hash = Column(String)

    pd = Column(SerializedDiagram)

    def __init__(self, pdcode):
        self.uid = pdcode.uid
        self.n_cross = pdcode.ncross
        self.n_comps = pdcode.ncomps
        self.hash = pdcode.hash
        self.pd = pdcode

    def __repr__(self):
        return "<Shadow(n_cross=%s, n_comps=%s, uid=%s, hash='%s')>"%(
            self.n_cross, self.n_comps, self.uid, self.hash)

class DiagramClass(Base):
    __tablename__ = 'diagram_classes'

    id = Column(Integer, primary_key=True)
    shadow_id = Column(Integer, ForeignKey('shadows.id'), index=True)
    shadow = relationship("Shadow", backref=backref('diagram_classes', order_by=id))

    homfly = Column(HOMFLYType, index=True)
    linktype_id = Column(Integer, ForeignKey('factorizations.id'), index=True)
    linktype = relationship("LinkFactorization", backref=backref('diagram_classes', order_by=id))

    diagrams = relationship("Diagram", backref="iso_class")

    def __repr__(self):
        return "<DiagramClass(n_cross=%s, n_comps=%s, shadow=%s, id=%s)>"%(
            self.shadow.n_cross, self.shadow.n_comps, self.shadow.uid, self.id)

class Diagram(Base):
    __tablename__ = 'diagrams'

    id = Column(Integer, primary_key=True)
    iso_class_id = Column(Integer, ForeignKey('diagram_classes.id'), index=True)

    cross_mask = Column(Mask)
    comp_mask = Column(Mask)

    @property
    def cross_mask_list(self):
        return int_to_bin_list(self.cross_mask, self.shadow.n_cross)
    @property
    def comp_mask_list(self):
        return int_to_bin_list(self.comp_mask, self.shadow.n_comps)

    @property
    def shadow(self):
        return self.iso_class.shadow

    @property
    def pd(self):
        pd = self.shadow.pd.copy()
        unsignature = (1-x for x in self.comp_mask_list)
        for comp in compress(range(pd.ncomps), unsignature):
           pd.reorient_component(comp, 0)
        pd.set_all_crossing_signs(self.cross_mask_list)
        return pd

    def __repr__(self):
        return "<Diagram(n_cross=%s, n_comps=%s, shadow=%s, cross_mask=%s, comp_mask=%s)>"%(
            self.shadow.n_cross, self.shadow.n_comps, self.shadow.uid,
            self.cross_mask, self.comp_mask)

from itertools import groupby
def dumb_get_factorization(session, facstr):
    """ facstr is of fmt "Knot[3, 1]#Knot[3, 1]#Knot[3, 1]*" """
    fzn = [(fact, len(list(grouper))) for
           fact,grouper in groupby(sorted(facstr.split("#")))]
    # fzn = ["Knot[3, 1]", 2, ...]

    q = session.query(LinkFactorization)

    # Use a join w/ correlated subquery to count num factors
    cfq = session.query(FactorizationFactor.product_id,
                        func.count(FactorizationFactor.product_id).label('fzn_count')).\
                        group_by(FactorizationFactor.product_id).subquery()
    q = q.join(cfq, cfq.c.product_id == LinkFactorization.id).filter(cfq.c.fzn_count==len(fzn))

    for factor, count in fzn:
        mirrored = False
        if "*" in factor:
            factor = factor[:-1]
            mirrored = True

        # Filter by this factor
        q = q.filter(LinkFactorization.factors.any(and_(
            FactorizationFactor.factor_id == session.query(LinkFactor.id)\
            .filter(LinkFactor.name == factor).filter(LinkFactor.mirrored == mirrored),
            FactorizationFactor.multiplicity == count
        )))

    # It is an error if there is more than 1 of this LF
    return q.one()

from sqlalchemy.ext.hybrid import Comparator
class FactorizationNameTransformer(Comparator):
    def operate(self, op, other):
        def transform(q):
            fzn = [(fact, len(list(grouper))) for
                   fact,grouper in groupby(sorted(other.split("#")))]

            cls = self.__clause_element__()

            cfq = q.session.query(
                FactorizationFactor.product_id,
                func.count(FactorizationFactor.product_id).label('fzn_count')).\
                group_by(FactorizationFactor.product_id).subquery()
            q = q.join(cfq, cfq.c.product_id == cls.id)\
                 .filter(cfq.c.fzn_count==len(fzn))

            for factor, count in fzn:
                mirrored = False
                if "*" in factor:
                    factor = factor[:-1]
                    mirrored = True

                # Filter by this factor
                q = q.filter(cls.factors.any(and_(
                    FactorizationFactor.factor_id == q.session.query(LinkFactor.id)\
                    .filter(LinkFactor.name == factor).filter(LinkFactor.mirrored == mirrored),
                    FactorizationFactor.multiplicity == count
                )))
            return q
        return transform

class LinkFactorization(Base):
    __tablename__ = 'factorizations'

    id = Column(Integer, primary_key=True)
    homfly = Column(HOMFLYType)

    factors = relationship('FactorizationFactor')

    n_splits = Column(Integer)
    n_cross = Column(Integer, index=True)
    n_comps = Column(Integer, index=True)

    @hybrid_property
    def name(self):
        return "#".join(
            "#".join(str(fassoc.factor) for _ in range(fassoc.multiplicity))
            for fassoc in self.factor)

    @name.comparator
    def name(cls):
        return FactorizationNameTransformer(cls)

    def __str__(self):
        return "#".join(
            "#".join(str(fassoc.factor) for _ in range(fassoc.multiplicity))
            for fassoc in self.factors)

    def __repr__(self):
        return "#".join(
            "#".join(repr(fassoc.factor) for _ in range(fassoc.multiplicity))
            for fassoc in self.factors)

class LinkFactor(Base):
    __tablename__ = 'factors'

    id = Column(Integer, primary_key=True)
    homfly = Column(HOMFLYType)
    name = Column(String)

    n_cross = Column(Integer)
    n_comps = Column(Integer)

    pd = Column(SerializedDiagram)

    mirrored = Column(Boolean)
    comp_mask = Column(Mask)

    @property
    def comp_mask_list(self):
        return int_to_bin_list(self.comp_mask, self.n_comps)

    def __str__(self):
        return "%s%s%s"%(
            self.name,
            "*" if self.mirrored else "",
            "_"+"".join(str(j) for j in self.comp_mask_list[1:]) if self.n_comps > 1 else "")


    def __repr__(self):
        return "<LinkFactor(%s%s%s)>"%(
            self.name,
            "*" if self.mirrored else "",
            "_"+"".join(str(j) for j in self.comp_mask_list[1:]) if self.n_comps > 1 else "")

#Diagram.factorizations = column_property(select([LinkFactorization])\
#                                         .where(LinkFactorization.n_cross==Shadow.n_cross)\
#                                         .where(LinkFactorization.homfly==Diagram.homfly))

def classify(session, diagram):
    return session.query(LinkFactorization)\
                  .filter(LinkFactorization.n_cross <= diagram.shadow.n_cross)\
                  .filter(LinkFactorization.homfly == diagram.homfly).all()

Base.metadata.create_all(engine)
Session = sessionmaker(bind=engine)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Load PD Shadows into the database.")
    parser.add_argument('shadows', metavar='N', type=int, nargs='+',
                        help="crossing counts for files to load in")
    args = parser.parse_args()

    session = Session()

    from libpl.pdstor import *

    for N in args.shadows:
        print "Adding shadows from ../../data/pdstors/%s.pdstor to database..."%N
        with open("../../data/pdstors/%s.pdstor"%N, "rb") as f:
            in_db = set((shadow.uid for shadow in (session.query(Shadow.uid)
                                                   .filter(Shadow.n_cross==N))))
            for shadow in PlanarDiagram.read_all(f, read_header=True):
                if shadow.ncomps > 1:
                    continue
                #(ret, ), = session.query(exists().where(
                #    Shadow.uid==shadow.uid).where(
                #        Shadow.n_cross==shadow.ncross))
                if shadow.uid in in_db:
                    continue
                db_shadow = Shadow(
                    pdcode = shadow
                )
                session.add(db_shadow)

    session.commit()
    session.close()
