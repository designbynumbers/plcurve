from libpl.pdcode import PlanarDiagram
from libpl import pdcode
from spherogram.links.links import Link, CrossingStrand, Strand, Crossing, CyclicList, CrossingEntryPoint
from spherogram.links.orthogonal import OrthogonalLinkDiagram, Face

import networkx, random, string

class PlanarCrossingStrand(CrossingStrand):
    def __init__(self, crossing, entry_point):
        super(PlanarCrossingStrand, self).__init__(crossing, entry_point)
        if isinstance(crossing, Strand):
            self.edge = crossing.label[0]
        else:
            zero_pos, _ = crossing.label.understrand_pos()
            self.edge = crossing.label[(entry_point+zero_pos)%4]
    def rotate(self, s=1):
        return PlanarCrossingStrand(self.crossing,
                                    (self.entry_point +s)%len(self.crossing.adjacent))
    def opposite(self):
        return PlanarCrossingStrand(*self.crossing.adjacent[self.entry_point])

class PlanarEntryPoint(PlanarCrossingStrand, CrossingEntryPoint):
    def next(self):
        c, e = self.crossing, self.entry_point
        s = 1 if isinstance(c, Strand) else 2
        return PlanarEntryPoint(*c.adjacent[ (e + s) % (2*s)])

    def other(self):
        nonzero_entry_point = 1 if self.crossing.sign == -1 else 3
        other = nonzero_entry_point if self.entry_point == 0 else 0
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
    strands[0][0] = tail.crossing[tail.entry_point]
    for i in range(n - 1):
        strands[i][1] = strands[i+1][0]
    strands[-1][1] = head.crossing[head.entry_point]

class OrthogonalPlanarDiagram(OrthogonalLinkDiagram):
    def __init__(self, pd):
        self.pd = pd = pd.copy()
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

    def ascii_data(self):
        """
        Returns:
        * a list of vertex positions
        * a list of arrows joining vertices
        * a list of crossings in the format (arrow over, arrow under)
        """
        emb = self.orthogonal_rep().basic_grid_embedding()
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
            vertex_positions.append( (10*(a+1), 10*(b+1)) )

        vert_indices = dict( (v,i) for i, v in enumerate(self.strand_CEPs))
        arrows, crossings = self.break_into_arrows()
        arrows = [ (vert_indices[a[0]], vert_indices[a[-1]], a) for a in arrows]

        return vertex_positions, arrows, crossings

### Orthogonal stuff ends here, start ascii stuff

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

class Pixel(object):
    def __init__(self, pixel, arrow=None):
        self.pixel = pixel
        self.arrow = arrow
    def __str__(self):
        return self.pixel
    def __repr__(self):
        return self.pixel

class Grid(list):
    def __init__(self, x, y):
        return super(Grid,self).__init__([[Pixel(" ")]*x for _ in range(y)])#[]
    def __str__(self):
        return "\n".join("".join(str(x) for x in row) for row in self)
    def __repr__(self):
        return "\n".join("".join(str(x) for x in row) for row in self)

def new_grid(x, y):
    return Grid(x,y)

def draw_label(grid,label,x,y,arrow=None,anchor="right"):
    if anchor == "left":
        oset = 0
    elif anchor == "center":
        oset = len(label)//2
    else:
        oset = len(label)-1
    for i,char in enumerate(label):
        crs_label_pixel = Pixel(char, arrow)
        grid[y][x-oset+i] = crs_label_pixel


def draw(orth):
    XBUFF = 3
    GRIDH = 4
    GRIDW = 14
    #orth = OrthogonalLinkDiagram(Link(K.pdcode()))
    #orth = OrthogonalLinkDiagram(K.as_spherogram())
    o_arrows, o_xings = orth.break_into_arrows()
    verts, arrows, xings = orth.ascii_data()
    verts = [((x/10-1)*GRIDW+XBUFF,(y/10-1)*GRIDH) for x,y in verts]
    arrow_pos = [(verts[arrow[0]], verts[arrow[1]], arrow[2]) for arrow in arrows]
    xings = [(o,u) for o,u,_ in xings]

    crs_pixels = {}
    edge_pixels = {}
    face_pixels = {}

    real_arrows = [arrow[2] for arrow in arrows]

    #print xings
    #print "----------"
    #for i, arrow in enumerate(arrow_pos):
    #    print i, arrow
    #print "~~~~~~~~~~"

    xmax = max([x for x,_ in verts])
    ymax = max([y for _,y in verts])
    grid = new_grid(xmax+1, ymax+1)
    for i, arrow in enumerate(arrow_pos):
        real_arrow = arrow[2]
        arrow = arrow[:2]
        points, head, pixel, headpixel, asign = arrow_to_pixels(*arrow)
        pixel = Pixel(pixel, real_arrow)
        headpixel = Pixel(headpixel, real_arrow)

        #if i not in (1,6):
        #    continue
        #print real_arrow

        k = 1
        for x,y in points:
            try:
                j = real_arrows.index(grid[y][x].arrow)
            except:
                j = None
            if j is None:
                grid[y][x] = pixel
            else:
                #a_xing = k if asign else -k-1
                #print o_arrows[i][a_xing].crossing, o_arrows[i][a_xing].crossing.adjacent
                #crs_pixels[real_arrow[a_xing].crossing] = (x-1,y-1)
                real_xings = [xing.crossing.label.index for xing in real_arrow[1:-1]]
                other_xings = [xing.crossing.label.index for xing in arrows[j][2][1:-1]]
                common_xing, = [xing for xing in real_xings if xing in other_xings]
                xing_strand_a = real_arrow[real_xings.index(common_xing)]
                #print xing_strand_a
                crs_pixels[common_xing] = (x-1, y-1)
                if (i, j) in xings:
                    #crs_pixels[xings.index((i,j))] = (x-1,y-1)
                    grid[y][x] = pixel
                else:
                    #crs_pixels[xings.index((j,i))] = (x-1,y-1)
                    assert (j, i) in xings
                k+=1
        x,y = head
        assert grid[y][x].arrow == real_arrow
        grid[y][x] = headpixel

    for n, (x,y) in crs_pixels.iteritems():
        draw_label(grid, "x"+str(n), x+2,y+2, n, anchor="left")

    for n, (xx,xy) in crs_pixels.iteritems():
        hz = grid[xy+1][xx+5].arrow
        vt = grid[xy-2][xx+1].arrow
        if hz:
            hz_xings = [xing.crossing.label.index for xing in hz]
            strand = hz[hz_xings.index(n)]
            #print "Strand", strand.edge
            hz_i = real_arrows.index(hz)
            #print arrow_pos[hz_i]
            ltr = arrow_pos[hz_i][0][0] < arrow_pos[hz_i][1][0]
            #label = str(strand.edge)+"--"+str(strand.edge.index)
            label = str(strand.edge.index)
            if not ltr:
                draw_label(grid, label, xx+5, xy+1, anchor="left")
            else:
                draw_label(grid, label, xx-3, xy+1, anchor="right")
        if vt:
            vt_xings = [xing.crossing.label.index for xing in vt]
            strand = vt[vt_xings.index(n)]
            #print "Strand", strand.edge
            vt_i = real_arrows.index(vt)
            #print arrow_pos[vt_i]
            ttb = arrow_pos[vt_i][0][1] < arrow_pos[vt_i][1][1]
            #label = str(strand.edge)+"  "+str(strand.edge.index)
            label = str(strand.edge.index)
            if not ttb:
                draw_label(grid, label, xx+1, xy+3, anchor="center")
            else:
                draw_label(grid, label, xx+1, xy-1, anchor="center")

    # Faces
    face_pos = dict()
    for n, (x, y) in crs_pixels.iteritems():
        strands = [PlanarCrossingStrand(orth.crossings[n], i) for i in range(4)]
        dirs = [orth.orientations[strand] for strand in strands]

        checks = [('right', (4,-1)), ('up', (-4, -1)), ('left', (-4,3)), ('down', (4, 3))]
        for dir, (dx, dy) in checks:
            strand = strands[dirs.index(dir)]
            face, = [i for i,F in enumerate(orth) if strand in F]
            face_pos[face] = (x+dx, y+dy)

    for i, (x,y) in face_pos.iteritems():
        draw_label(grid, "f"+str(i), x, y, anchor="left")

    print grid#"\n".join("".join(str(x) for x in row) for row in grid)
    #print K.edges
    #print xings

    #import IPython
    #IPython.embed()

if __name__ == "__main__":
    #K = PlanarDiagram.simple_chain(3)
    K = PlanarDiagram.torus_knot(2,8)
    print repr(K)
    print K.faces
    draw(OrthogonalPlanarDiagram(K))
