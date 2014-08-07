from libpl.pdstor import PDDatabase
from cassandra.cluster import Cluster
import cPickle

DEFAULT_PATH="../../data/pdstors"
class CassandraDatabase(PDDatabase):
    def __init__(self, crossings_list=[3,4], dirloc=DEFAULT_PATH, debug=False):
        self.cluster = Cluster(['23.23.133.228'])
        self.session = self.cluster.connect("diagrams")
        self.queries = {}
        for ncross in crossings_list:
            self.queries[ncross] = self.session.prepare("""INSERT INTO cross%s
            (pdcode, homfly, pythonobject, topohash) VALUES
            (?, ?, ?, ?);"""%ncross)

        super(CassandraDatabase, self).__init__(
            crossings_list, dirloc,
            amortize=True, debug=debug,
            callbacks=[self.insert,])

    def insert(self, pd, homfly, uid):
        print homfly
        self.futures.append(
            self.session.execute_async(self.queries[pd.ncross],
                                       (str(uid), homfly, cPickle.dumps(pd,2), pd.hash))
        )
