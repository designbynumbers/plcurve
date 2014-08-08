from libpl.pdstor import PDDatabase
from cassandra.cluster import Cluster
import cPickle

DEFAULT_PATH="../../data/pdstors"
class CassandraDatabase(PDStoreExpander):
    def __init__(self, hosts=['localhost'], dirloc=DEFAULT_PATH, debug=False):
        if isinstance(hosts, basestring):
            hosts = [hosts]
        self.cluster = Cluster(hosts)
        self.session = self.cluster.connect("diagrams")
        self.queries = {}
        self.futures = []
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
