from cassandra.cluster import Cluster

COLD_IP = "23.23.133.228"

class ColdInterface(object):
    def __init__(self, host):
        self._cluster = Cluster(host)
        self._session = None

    @property
    def session(self):
        if not self._session:
            self.connect("diagrams")
        return self._session

    def _prepare(self, stmt):
        if stmt not in self._statements:
            self._statements[stmt] = session.prepare(stmt)

        return self._statements[stmt]

    def connect(self):
        self._session = self._cluster.connect()
        self._statements = dict()

    def pd_by_id(self, uid):
        ncross, filepos, crossmask, compmask = uid
        stmt = self._prepare("SELECT * FROM cross%s WHERE pdcode=?" % ncross)
        pds = self.session.execute(stmt, [str(uid)])
