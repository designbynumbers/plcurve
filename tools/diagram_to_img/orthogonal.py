"""
orthogonal.py

Provides the OrthogonalPlanarDiagram class, which represents an
orthogonal representation of a PlanarDiagram object on an XY-plane.

We get pretty dirty in `spherogram` class inheritance here, and we
have to be careful as `spherogram` changes versions. Currently, this
works fine with `spherogram` version 1.6.1 (hchapman 07/05/2017)
"""

from libpl.pdcode import PlanarDiagram
from libpl import pdcode
from spherogram.links.links import Link, CrossingStrand, Strand, Crossing, CrossingEntryPoint
from spherogram.links.orthogonal import OrthogonalLinkDiagram, Face

import networkx, random, string

class PlanarCrossingStrand(CrossingStrand):
    def __init__(self, crossing, strand_index):
        super(PlanarCrossingStrand, self).__init__(crossing, strand_index)
        if isinstance(crossing, Strand):
            self.edge = crossing.label[0]
        else:
            zero_pos, _ = crossing.label.understrand_pos()
            self.edge = crossing.label[(strand_index+zero_pos)%4]
    def rotate(self, s=1):
        return PlanarCrossingStrand(self.crossing,
                                    (self.strand_index +s)%len(self.crossing.adjacent))
    def opposite(self):
        return PlanarCrossingStrand(*self.crossing.adjacent[self.strand_index])

class PlanarEntryPoint(PlanarCrossingStrand, CrossingEntryPoint):
    def next(self):
        c, e = self.crossing, self.strand_index
        s = 1 if isinstance(c, Strand) else 2
        return PlanarEntryPoint(*c.adjacent[ (e + s) % (2*s)])

    def other(self):
        nonzero_strand_index = 1 if self.crossing.sign == -1 else 3
        other = nonzero_strand_index if self.strand_index == 0 else 0
        return PlanarEntryPoint(self.crossing, other)

class PlanarFace(Face):
    def __init__(self, crossings, face, exterior=False):
        self.raw_edges = list(reversed([edge for edge in zip(face, face.signs)]))
        heads = [(edge.head, edge.headpos) if sign == 0 else (edge.tail, edge.tailpos) for
                 edge,sign in self.raw_edges]
        strands = [PlanarCrossingStrand(crossings[i], (4+pos-crossings[i].label.understrand_pos()[0])%4)
                   for (i,pos),edge in zip(heads, self.raw_edges)]
        super(PlanarFace, self).__init__(None, strands)

    def make_exterior(self):
        new_strands = [strand.opposite() for strand in self]
        super(PlanarFace, self).__init__(None, list(new_strands), exterior=True)

class PdCrossing(Crossing):
    def __str__(self):
        return "PD %s%s"%(self.label.index, "+" if self.label.sign else "-")
    def __repr__(self):
        return str(self)

class PlanarStrand(Strand):
    def __repr__(self):
        return "%s %s"%self.label

def subdivide_edge(crossing_strand, n):
    """
    Given a PlanarCrossingStrand, subdivides the edge for which it's the *head*
    into (n + 1) pieces.

    WARNING: this breaks several of the link's internal data structures.
    """
    head = crossing_strand
    backwards = not (head in head.crossing.entry_points())
    if backwards:
        head = head.opposite()
    tail = head.opposite()

    strands = [PlanarStrand(label=(head.edge, i)) for i in range(n)]
    strands[0][0] = tail.crossing[tail.strand_index]
    for i in range(n - 1):
        strands[i][1] = strands[i+1][0]
    strands[-1][1] = head.crossing[head.strand_index]

class OrthogonalPlanarDiagram(OrthogonalLinkDiagram):
    def __init__(self, pd):
        #pd = pd.copy()
        self.pd = pd
        self.crossings = [PdCrossing(label=x) for x in pd.crossings]
        for x, cross in zip(pd.crossings, self.crossings):
            # Identify which of our PD indices is SGM index 0:
            zero_index, _ = x.understrand_pos()
            for pd_pos in range(4):
                adj_x, adj_pd_pos = x.adjacent(pd_pos)
                adj_zero_index, _ = adj_x.understrand_pos()
                i = (4+pd_pos-zero_index)%4
                adj_pos = (4+adj_pd_pos-adj_zero_index)%4
                cross.adjacent[i] = self.crossings[adj_x.index][adj_pos]
            cross.sign = x.sign*2-1

        list.__init__(self, [PlanarFace(self.crossings, F) for F in pd.faces])
        exterior_face = max(self, key=len)
        exterior_face.exterior = True
        self.face_network = self.flow_networkx()
        self.bend()
        self.orient_edges()
        self.edges = sum([F for F in self], [])
        self.repair_components()

    def bend(self):
        """
        Computes a minimal size set of edge bends that allows the link diagram
        to be embedded orthogonally. This follows directly Tamassia's first
        paper.
        """
        N = self.face_network
        flow = networkx.min_cost_flow(N)
        for a, flows in flow.iteritems():
            for b, w_a in flows.iteritems():
                if w_a and set(['s', 't']).isdisjoint(set([a, b])):
                    w_b = flow[b][a]
                    A, B = self[a], self[b]
                    e_a, e_b = A.edge_of_intersection(B)
                    turns_a = w_a*[1] + w_b*[-1]
                    turns_b = w_b*[1] + w_a*[-1]
                    subdivide_edge(e_a, len(turns_a))
                    A.bend(e_a, turns_a)
                    B.bend(e_b, turns_b)

    def repair_components(self):
        """???"""
        components = [PlanarEntryPoint(
            self.crossings[component[0].next_crossing().index],
            (4+component[0].headpos-component[0].next_crossing().understrand_pos()[0])%4).component()
                      for component in self.pd.components]
        self.strand_CEP_to_component = stc = dict()
        self.strand_CEPs = []
        for n, component in enumerate(components):
            for ce in component:
                if isinstance(ce.crossing, Strand):
                    stc[ce] = n
                    self.strand_CEPs.append(ce)

    def ascii_data(self, emb, xscale=1, yscale=1):
        """
        Returns:
        * a list of vertex positions
        * a list of arrows joining vertices
        * a list of crossings in the format (arrow over, arrow under)
        """
        x_max = max(a for a,b in emb.values())
        y_max = max(b for a,b in emb.values())

        # We rotate things so the long direction is horizontal.  The
        # Plink canvas coordinate system forces us to flip things to
        # preserve crossing type.
        vertex_positions = []
        for v in self.strand_CEPs:
            #if x_max >= y_max:
            a, b = emb[v.crossing]
            b = y_max - b
            #else:
            #b, a = emb[v.crossing]
            vertex_positions.append( (xscale*(a+1), yscale*(b+1)) )

        vert_indices = dict( (v,i) for i, v in enumerate(self.strand_CEPs))
        arrows, crossings = self.break_into_arrows()
        arrows = [ (vert_indices[a[0]], vert_indices[a[-1]], a) for a in arrows]

        return vertex_positions, arrows, crossings

def arrow_to_pixels(tail, head):#head, tail):
    horizontal = bool(tail[0] - head[0])
    z_i = 0 if horizontal else 1
    z = range(1+min(tail[z_i], head[z_i]), max(tail[z_i], head[z_i]))
    if horizontal:
        y = tail[1]
        ltr = head[0]>tail[0]
        return (((x,y) for x in z),
                (head[0]-1 if ltr else head[0]+1, head[1]),
                "-", ">" if ltr else "<",
                not ltr)
    else:
        x = tail[0]
        ttb = head[1]>tail[1]
        return (((x,y) for y in z),
                (head[0], head[1]-1 if ttb else head[1]+1),
                "|", "v" if ttb else "^",
                ttb)
