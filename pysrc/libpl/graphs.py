SIGN_CHAR = {0: "-", 1: "+", 2: "?"}

class DiEdge(object):
    def __init__(self, tail, head):
        self.tail = tail
        self.head = head

    def __repr__(self):
        return "DiEdge(%s,%s)"%(self.tail, self.head)
    def __str__(self):
        return "%s->%s"%(self.tail, self.head)

class Face(object):
    def __init__(self, cycle):
        self.cycle = cycle

    def __repr__(self):
        return "Face(%s)"%repr(self.cycle)
    def __str__(self):
        return str(self.cycle)

class SignedFace(Face):
    def __init__(self, cycle, sign):
        super(SignedFace, self).__init__(cycle)
        self.sign = sign
        
    def __repr__(self):
        return "SignedFace(%s,%s)"%(repr(self.cycle), self.sign)
    def __str__(self):
        return "%s%s"%(" ".join(str(e) for e in self.cycle), SIGN_CHAR[self.sign])

class PlanarSignedFaceDigraph(object):
    def __init__(self):
        self.verts = []
        self.edges = []
        self.faces = []

    def add_edge(self, tail, head):
        self.edges.append(DiEdge(tail, head))

    def add_face(self, edge_cycle, sign):
        self.faces.append(SignedFace(edge_cycle, sign))