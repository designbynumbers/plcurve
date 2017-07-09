from libpl.pdcode import PlanarDiagram
from libpl import pdcode
from spherogram.links.links import Link, CrossingStrand, Strand, Crossing, CrossingEntryPoint
from spherogram.links.orthogonal import OrthogonalLinkDiagram, Face

import cairo
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

    def ascii_data(self, xscale=1, yscale=1):
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
            vertex_positions.append( (10*xscale*(a+1), 10*yscale*(b+1)) )

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
        super(Grid,self).__init__([[Pixel(" ")]*x for _ in range(y)])
    def __str__(self):
        return "\n".join("".join(str(x) for x in row) for row in self)
    def __repr__(self):
        return "\n".join("".join(str(x) for x in row) for row in self)

class Grid2(list):
    def __init__(self, x, y, w=3, h=3):
        self._w, self._h = w+1,h+1
        super(Grid2, self).__init__([[Pixel(" ")]*(self._w*x) for _ in range(self._h*y)])

    def add_vert(self, point, char):
        x,y = point
        self[self._h*y][self._w*x] = Pixel(char)

    def plot_edge(self, tail, head):
        x1,y1 = tail
        x2,y2 = head
        dx,dy = x1-x2, y1-y2
        #print tail,head, dx, dy
        if (dx and dy) or (not dx and not dy):
            raise Exception("Can only plot horizontal or vertical edges")

        if dx: # Horizontal edge
            y = self._h*y1
            if dx < 0: # Edge is left-to-right
                for x in range(self._w*x1+1, self._w*x2-1):
                    self[y][x] = Pixel("-")
                self[y][self._w*x2-1] = Pixel(">")
            else:
                for x in range(self._w*x2+1, self._w*x1-1):
                    self[y][x] = Pixel("-")
                self[y][self._w*x1-1] = Pixel(">")

    def display(self, h=1, w=1):
        print "\n".join("".join(str(x) for x in row) for row in self)

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

def draw2(orth):
    orep = orth.orthogonal_rep()
    emb = orep.basic_grid_embedding()

    grid = Grid2(5,5)

    for crossing, point in emb.iteritems():
        grid.add_vert(point, "x")

    for tail_x, head_x in orep.edges:
        tail_p, head_p = emb[tail_x], emb[head_x]
        grid.plot_edge(tail_p, head_p)

    grid.display()

def draw(orth, skip_start=2, skip_end=4):
    XBUFF = 3
    GRIDH = 4
    GRIDW = 14
    #orth = OrthogonalLinkDiagram(Link(K.pdcode()))
    #orth = OrthogonalLinkDiagram(K.as_spherogram())
    o_arrows, o_xings = orth.break_into_arrows()
    verts, arrows, xings = orth.ascii_data()
    verts = [((x/10-1)*GRIDW+XBUFF,(y/10-1)*GRIDH) for x,y in verts]
    arrow_pos = [(verts[arrow[0]], verts[arrow[1]], arrow[2]) for arrow in arrows]
    #print xings
    xings = [(o,u) for o,u,_,_ in xings]

    crs_pixels = {}
    edge_pixels = {}
    face_pixels = {}

    real_arrows = [arrow[2] for arrow in arrows]

    ##print xings
    ##print "----------"
    #for i, arrow in enumerate(arrow_pos):
    #    #print i, arrow
    ##print "~~~~~~~~~~"

    xmax = max([x for x,_ in verts])
    ymax = max([y for _,y in verts])
    grid = new_grid(xmax+1, ymax+1)
    for i, arrow in enumerate(arrow_pos):
        ##print arrow
        real_arrow = arrow[2]
        arrow = arrow[:2]
        points, head, pixel, headpixel, asign = arrow_to_pixels(*arrow)
        pixel = Pixel(pixel, real_arrow)
        headpixel = Pixel(headpixel, real_arrow)

        if len(real_arrow) == 2 and skip_start <= real_arrow[0].edge.index <= skip_end:
            # If real_arrow is len 2, there's no interesting behavior
            # so if we want to skip this guy, we certainly can
            continue

        k = 1
        for x,y in points:
            try:
                # Is there already an arrow at this pixel?
                j = real_arrows.index(grid[y][x].arrow)
            except:
                j = None

            if j is None:
                # There is no arrow at this pixel! We have free reign

                grid[y][x] = pixel

            else:
                # There is already an arrow at this pixel.
                # Hence we must be crossing over or under.
                real_xings = [xing.crossing.label.index for xing in real_arrow[1:-1]]
                other_xings = [xing.crossing.label.index for xing in arrows[j][2][1:-1]]
                common_xing, = [xing for xing in real_xings if xing in other_xings]
                xing_strand_a = real_arrow[real_xings.index(common_xing)]
                #print pixel, xing_strand_a.edge.index
                #pixel = Pixel(str(xing_strand_a.edge.index), real_arrow)

                ##print xing_strand_a
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

            hz_i = real_arrows.index(hz)

            # Is arrow left-to-right?
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
            ##print "Strand", strand.edge
            vt_i = real_arrows.index(vt)
            ##print arrow_pos[vt_i]
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

    #print grid#"\n".join("".join(str(x) for x in row) for row in grid)
    ##print K.edges
    ##print xings

    #import IPython
    #IPython.embed()

from collections import defaultdict

def intersect_arrows(over, under, verts):
    o_si, o_ei = over
    u_si, u_ei = under
    over_start, over_end = verts[o_si], verts[o_ei]
    under_start, under_end = verts[u_si], verts[u_ei]
    #print "~~~~"
    #print over_start, over_end, under_start, under_end
    if over_start[0] == over_end[0]:
        # over is horizontal
        x = over_start[0]
        y = under_start[1]
    else:
        x = under_start[0]
        y = over_start[1]
    return (x, y)

def draw_orthogonal_projection(orth, fname,
                               width=100., height=100.,
                               skip_start=0, skip_end=0,
                               cross_width=4):
    XBUFF = 3
    GRIDH = 5
    GRIDW = 5

    ## Natural scale is a grid by 10s
    xscale = 2.
    yscale = 2.

    verts, arrows, xings = orth.ascii_data(1., 1.)

    maxx = max(v[0] for v in verts)
    maxy = max(v[1] for v in verts)

    sqwidth = min(width, height)
    ctr_x, ctr_y = width/2, height/2
    rad = sqwidth/2-30

    scale = min(maxx, maxy)
    xscale = width/(maxx+10)
    yscale = height/(maxy+10)
    ssc = min(xscale, yscale)

    verts, arrows, xings = orth.ascii_data(ssc, ssc)


    # Goal 1 is to break up arrows into crossing chunks
    # *------|------|------|------>* should become...
    # *------|--)(--|--)(--|------>*
    # so that an arrow piece either has
    #   a) NO crossings, and is part of precisely one edge, or
    #   b) 1 crossing, and is part of two different edges
    # for planar diagrams, edges in case b) are guaranteed to be distinct

    easy_arrows = []
    broken_arrows = []

    cross_map = defaultdict(list)
    for over, under, _, crs in xings:
        #print _, crs
        #print arrows[over]
        x,y = intersect_arrows(arrows[over][:2], arrows[under][:2], verts)
        #print x, y

        # cross map~ arrow: (is_over, otherarrow, intersection, crs)
        cross_map[over].append((True, under, (x, y), crs))
        cross_map[under].append((False, over, (x, y), crs))

    for a_i, arrow in enumerate(arrows):
        #print arrow[2]
        if len(arrow[2]) == 2:
            # only corner pseudocrossings, no real crossings
            # so this arrow is an honest easy arrow
            easy_arrows.append((arrow[2][0].edge.index, verts[arrow[0]], verts[arrow[1]]))

        elif len(arrow[2]) == 3:
            # this arrow has precisely 1 crossing; it hence doesn't
            # need to be broken, but it does encompass two edges
            broken_arrows.append((arrow[2][1].edge.index, arrow[2][2].edge.index,
                                  cross_map[a_i][0][0], cross_map[a_i][0][2],
                                  verts[arrow[0]], verts[arrow[1]]))

        else:
            # this arrow has multiple crossings; we have to break it
            # at new pseudocrossings inbetween its crossings
            arrowtail = verts[arrow[0]]
            arrowhead = verts[arrow[1]]

            sgn = 1 if arrowtail[0]+arrowtail[1] < arrowhead[0]+arrowhead[1] else -1
            arrow_cross = sorted(cross_map[a_i], key=lambda t: sgn*(t[2][0]+t[2][1]))
            # crossings here are sorted from tail to head

            # now add in the breaks
            tail_vtx = arrowtail
            for i in range(len(arrow_cross)-1):
                # there is one fewer break than crossings
                # break vertex is midpoint of next two crossings
                #print arrow_cross
                head_vtx = ((arrow_cross[i][2][0]+arrow_cross[i+1][2][0])/2,
                            (arrow_cross[i][2][1]+arrow_cross[i+1][2][1])/2)
                tail_edge_i = arrow[2][i+1].edge.index
                head_edge_i = arrow[2][i+2].edge.index
                is_over, _, crspos, crs = arrow_cross[i]
                broken_arrows.append((tail_edge_i, head_edge_i,
                                      is_over, crspos, tail_vtx, head_vtx))
                #print "brk", self.broken_arrows
                tail_vtx = head_vtx
            head_vtx = arrowhead
            tail_edge_i = arrow[2][-2].edge.index
            head_edge_i = arrow[2][-1].edge.index
            is_over, _, crspos, crs = arrow_cross[-1]
            broken_arrows.append((tail_edge_i, head_edge_i,
                                  is_over, crspos, tail_vtx, head_vtx))
            #print "brk", self.broken_arrows

    surf = cairo.SVGSurface(fname, width, height)
    cr = cairo.Context(surf)
    for edge, arrowtail, arrowhead in easy_arrows:
        # skip an easy arrow if its edge is skippable
        if skip_start <= edge <= skip_end:
            continue
        elif (skip_end < skip_start and
              (edge <= skip_end or edge >= skip_start)):
            continue

        cr.move_to(*arrowtail)
        cr.line_to(*arrowhead)
        cr.stroke()

    for tailedge, headedge, is_over, crossloc, arrowtail, arrowhead in broken_arrows:
        # skip a broken arrow if both its tail and head are skippable
        if (skip_start <= tailedge <= skip_end and
            skip_start <= headedge <= skip_end):
            continue
        elif (skip_end < skip_start and (tailedge <= skip_end or skip_start <= tailedge) and
              (headedge <= skip_end or skip_start <= headedge)):
            continue

        #print "~~~", arrowtail, arrowhead
        is_horizontal = arrowtail[0] != arrowhead[0]

        if not is_over:
            if is_horizontal:
                sgn = 1 if arrowtail[0] < arrowhead[0] else -1
                #print "hzunder", arrowtail, crossloc, arrowhead
                cr.move_to(arrowtail[0], arrowtail[1])
                cr.line_to(crossloc[0] - sgn*cross_width, crossloc[1])
                cr.move_to(crossloc[0] + sgn*cross_width, crossloc[1])
                cr.line_to(arrowhead[0], arrowhead[1])
                cr.stroke()
            else:
                sgn = 1 if arrowtail[1] < arrowhead[1] else -1
                #print "vtunder", arrowhead, arrowtail, crossloc
                cr.move_to(arrowtail[0], arrowtail[1])
                cr.line_to(crossloc[0], crossloc[1] - sgn*cross_width)
                cr.move_to(crossloc[0], crossloc[1] + sgn*cross_width)
                cr.line_to(arrowhead[0], arrowhead[1])
                cr.stroke()
        else:
            cr.move_to(*arrowtail)
            cr.line_to(*arrowhead)
            cr.stroke()
            cr.arc(crossloc[0], crossloc[1], cross_width/2, 0, 2*pi)
            cr.fill()

    surf.write_to_png("knot.png")

unkn_hf = {}

def fuzzy_cls(hf):
    hf = str(hf)
    if hf in hf_to_kt:
        return hf_to_kt[hf]
    if hf in unkn_hf:
        return unkn_hf[hf]
    unkn_hf[hf] = "?_{}".format(len(unkn_hf)+1)
    return unkn_hf[hf]

from math import atan, asin, acos, cos, sin

import sys

#w.create_rectangle(0, 0, 500, 500, fill="red")

from math import sin, cos, pi, tan, asin
from time import sleep

def circle_point(center, radius, angle):
    x, y = center
    return (x + radius*cos(angle), y + radius*sin(angle))

def display_pdcode(circ_imm):
    pd = circ_imm.to_PlanarDiagram()
    #pd.edit_copy()
    #pd.as_spherogram().view()
    ##print pd

if __name__ == "__main__":
    #K = PlanarDiagram.simple_chain(2)
    #K = PlanarDiagram.torus_knot(2,5)
    K = PlanarDiagram.random_diagram(100, 3, 0)
    #with open("5_1slip.pd", "r") as f:
    #    K = PlanarDiagram.read(f)
    ##print repr(K)

    pd = K
    N = pd.ncross

    orth = OrthogonalPlanarDiagram(K)
    draw_orthogonal_projection(orth, "knot2.svg",height=400, width=400)
