from cassandra.cluster import Cluster
from cPickle import loads, dumps
from libpl.pdstor import PDStoreExpander
from libpl.pdcode import PlanarDiagram

DEFAULT_PATH = "../../data/pdstors"
KEYSPACE = "diagrams"

INSERT_TEMPLATE_NEW_CROSSN = """INSERT INTO cross%s
(pdcode, homfly, pythonobject, topohash) VALUES
(?, ?, ?, ?);"""

class ColdInterface(PDStoreExpander):
    def __init__(self, hosts=['localhost'], dirloc=DEFAULT_PATH, debug=False):
        if isinstance(hosts, basestring):
            hosts = [hosts]
        self._cluster = Cluster(hosts)
        self._session = None
        super(ColdInterface, self).__init__(dirloc=dirloc,amortize=True,debug=debug)

    @property
    def session(self):
        if not self._session:
            self.connect()
        return self._session

    def _prepare(self, stmt):
        session = self.session
        if stmt not in self._statements:
            self._statements[stmt] = session.prepare(stmt)

        return self._statements[stmt]

    def connect(self):
        self._session = self._cluster.connect(KEYSPACE)
        self._statements = dict()

    def insert_all(self, ncrosses):
        """insert_all(ncrosses) -> list of futures

        Insert all diagrams determined by the shadows with N crossings
        for N in ``ncrosses`` into the Cassandra database table crossN.
        """
        futures = []
        session = self.session
        statements = {N: self._prepare(INSERT_TEMPLATE_NEW_CROSSN % N) for
                      N in ncrosses}

        for pdc, homfly, pdid in self.open(ncrosses):
            futures.append(
                session.execute_async(statements[pdc.ncross],
                                      (str(pdid), homfly, dumps(pdc,2), pdc.hash))
            )
        return futures

    def insert_all_and_wait(self, ncrosses):
        """insert_all_and_wait(ncrosses)

        Like ``insert_all()`` but wait for all the queries to finish executing."""
        for future in self.insert_all(ncrosses):
            future.result()

    def pd_by_id(self, uid):
        ncross, filepos, crossmask, compmask = uid
        stmt = self._prepare("SELECT * FROM cross%s WHERE pdcode=?" % ncross)
        pds = self.session.execute(stmt, [str(uid)])
        if not pds:
            return None
        assert(len(pds) == 1)
        return loads(pds[0].pythonobject)
