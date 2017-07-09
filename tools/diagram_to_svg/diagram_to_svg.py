#!/usr/bin/python
from libpl.pdcode import PlanarDiagram
from libpl import pdcode

from collections import defaultdict

import cairo
import json

from orthogonal import OrthogonalPlanarDiagram

import os.path
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "config.json")) as f:
    colors = json.load(f)["colors"]
colors = [tuple(i/255. for i in color) for color in colors]

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

def get_orthogonal_arrows(orth, width=100., height=100.,
                          scale=None, xscale=None, yscale=None):
    """
    Using an OrthogonalPlanarDiagram create a list of primitive "arrow" shapes
    for drawing.
    """

    XBUFF = 3
    GRIDH = 5
    GRIDW = 5

    # Fix a grid embedding
    emb = orth.orthogonal_rep().basic_grid_embedding()

    # Case on whether we should scale to w/h or were given a scale
    if scale is not None:
        verts, arrows, xings = orth.ascii_data(emb, scale, scale)
        maxx = max(v[0] for v in verts)
        maxy = max(v[1] for v in verts)
        width, height = scale+maxx, scale+maxy

    elif xscale is not None and yscale is not None:
        verts, arrows, xings = orth.ascii_data(emb, xscale, yscale)
        maxx = max(v[0] for v in verts)
        maxy = max(v[1] for v in verts)
        width, height = xscale+maxx, yscale+maxy
    else:
        verts, arrows, xings = orth.ascii_data(emb, 1., 1.)

        maxx = max(v[0] for v in verts)
        maxy = max(v[1] for v in verts)

        sqwidth = min(width, height)
        ctr_x, ctr_y = width/2, height/2
        rad = sqwidth/2-30

        scale = min(maxx, maxy)
        xscale = width/(maxx+1)
        yscale = height/(maxy+1)
        ssc = min(xscale, yscale)

        verts, arrows, xings = orth.ascii_data(emb, xscale, yscale)

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
        # The format of arrow is:
        # (vert1_i, vert2_i, [CEP1, CEP2, ???])
        #print arrow[2]
        if len(arrow[2]) == 2:
            # only corner pseudocrossings, no real crossings
            # so this arrow is an honest easy arrow
            easy_arrows.append((arrow[2][0].edge.index, verts[arrow[0]], verts[arrow[1]],
                                arrow[2][0].edge.component_index_pos()[0]))

        elif len(arrow[2]) == 3:
            # this arrow has precisely 1 crossing; it hence doesn't
            # need to be broken, but it does encompass two edges
            broken_arrows.append((arrow[2][1].edge.index, arrow[2][2].edge.index,
                                  cross_map[a_i][0][0], cross_map[a_i][0][2],
                                  verts[arrow[0]], verts[arrow[1]],
                                  arrow[2][0].edge.component_index_pos()[0]))

        else:
            # this arrow has multiple crossings; we have to break it
            # at new pseudocrossings inbetween its crossings
            arrowtail = verts[arrow[0]]
            arrowhead = verts[arrow[1]]
            comp = arrow[2][0].edge.component_index_pos()[0]

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
                                      is_over, crspos, tail_vtx, head_vtx, comp))
                #print "brk", self.broken_arrows
                tail_vtx = head_vtx
            head_vtx = arrowhead
            tail_edge_i = arrow[2][-2].edge.index
            head_edge_i = arrow[2][-1].edge.index
            is_over, _, crspos, crs = arrow_cross[-1]
            broken_arrows.append((tail_edge_i, head_edge_i,
                                  is_over, crspos, tail_vtx, head_vtx, comp))
            #print "brk", self.broken_arrows

    return width, height, easy_arrows, broken_arrows

class DrawingContext(object):
    pass

class CairoDrawingContext(DrawingContext):
    def __init__(self, fname):
        import cairo
        self.fname = fname

    def open(self, width, height):
        self.surf = cairo.SVGSurface(self.fname, width, height)
        self.cr = cairo.Context(self.surf)
        self.cr.set_line_width(1.)

    def set_thickness(self, thickness):
        self.cr.set_line_width(thickness)

    def set_color(self, r, g, b):
        self.cr.set_source_rgb(r, g, b)

    def draw_line(self, start, end):
        self.cr.move_to(*start)
        self.cr.line_to(*end)
        self.cr.stroke()

    def draw_circle(self, center, radius):
        self.cr.arc(center[0], center[1], radius, 0, 2*pi)
        self.cr.fill()

    def close(self):
        pass

class PILDrawingContext(DrawingContext):
    def __init__(self, fname):
        import PIL.Image, PIL.ImageDraw
        self.fname = fname
        self.color = (0,0,0)
        self.thickness = 1

    def open(self, width, height):
        self.image = PIL.Image.new("RGB", (width, height), (255,255,255))
        self.draw = PIL.ImageDraw.Draw(self.image)

    def set_thickness(self, thickness):
        self.thickness = thickness

    def set_color(self, r, g, b):
        self.color = (int(r*255), int(g*255), int(b*255))

    def draw_line(self, start, end):
        self.draw.line((start, end), self.color, width=self.thickness)

    def draw_circle(self, center, radius):
        self.draw.ellipse([(center[0]-radius, center[1]-radius),
                           (center[0]+radius, center[1]+radius)],
                          self.color)

    def close(self):
        with open(self.fname+".png", "w") as f:
            self.image.save(f, "PNG")

def _draw_pd(pd, dc, skip_start=0, skip_end=0,
             cross_width=4, cross_radius=0,
             line_width=1, **kwargs):
    orth = OrthogonalPlanarDiagram(pd)
    width, height, easy_arrows, broken_arrows = get_orthogonal_arrows(
        orth, **kwargs)

    dc.open(width, height)
    print width,height

    dc.set_thickness(line_width)

    ## "Easy" arrows were never broken
    for edge, arrowtail, arrowhead, comp in easy_arrows:
        # skip an easy arrow if its edge is skippable
        if skip_start == 0 and skip_end == 0:
            pass
        elif skip_start <= edge <= skip_end:
            continue
        elif (skip_end < skip_start and
              (edge <= skip_end or edge >= skip_start)):
            continue

        dc.set_color(*colors[comp%len(colors)])
        dc.draw_line(arrowtail, arrowhead)

    ## "Broken" arrows were long arrows broken apart by crossings
    for tailedge, headedge, is_over, crossloc, arrowtail, arrowhead, comp in broken_arrows:
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
                dc.set_color(*colors[comp%len(colors)])
                dc.draw_line((arrowtail[0], arrowtail[1]),
                             (crossloc[0] - sgn*cross_width, crossloc[1]))
                dc.draw_line((crossloc[0] + sgn*cross_width, crossloc[1]),
                             (arrowhead[0], arrowhead[1]))
            else:
                sgn = 1 if arrowtail[1] < arrowhead[1] else -1
                #print "vtunder", arrowhead, arrowtail, crossloc
                dc.set_color(*colors[comp%len(colors)])
                dc.draw_line((arrowtail[0], arrowtail[1]),
                             (crossloc[0], crossloc[1] - sgn*cross_width))
                dc.draw_line((crossloc[0], crossloc[1] + sgn*cross_width),
                             (arrowhead[0], arrowhead[1]))

        else:
            dc.set_color(*colors[comp%len(colors)])
            dc.draw_line(arrowtail, arrowhead)

            if cross_radius:
                dc.draw_circle(crossloc, cross_radius)

    dc.close()

def draw_pd_cairo(pd, fname, **kwargs):
    _draw_pd(pd, CairoDrawingContext(fname),
             **kwargs)

def draw_pd_pil(pd, fname, **kwargs):
    _draw_pd(pd, PILDrawingContext(fname),
             **kwargs)

from math import pi

def _args_to_draw_kwargs(args):
    return {"height": args.height,
            "width": args.width,
            "scale": args.scale,
            "xscale": args.xscale,
            "yscale": args.yscale,
            "line_width": args.line_width,
            "cross_width": args.cross_width,
            "cross_radius": args.cross_radius}

def pdstor_command(args):
    for pd in PlanarDiagram.read_all(args.pdstor,
                                     read_header=args.read_header):
        pd.randomly_assign_crossings()
        args.draw_pd(pd, "diagram_{}.svg".format(pd.uid),
                     **_args_to_draw_kwargs(args))

def pdcode_command(args):
    i = 0
    pd = PlanarDiagram.read_knot_theory(args.pdcode)
    while pd is not None:
        args.draw_pd(pd, "diagram_pdcode_{}.svg".format(i),
                     **_args_to_draw_kwargs(args))
        pd = PlanarDiagram.read_knot_theory(args.pdcode)
        i += 1

def random_command(args):
    pd = PlanarDiagram.random_diagram(args.size)
    args.draw_pd(pd, "diagram_random.svg".format(pd.uid),
                 **_args_to_draw_kwargs(args))


if __name__ == "__main__":
    import sys, argparse

    parser = argparse.ArgumentParser(
        description='Generate SVG images from PlanarDiagrams')

    fmt_group = parser.add_mutually_exclusive_group()
    fmt_group.add_argument("--cairo-svg", action='store_const',
                           dest="draw_pd",
                           const=draw_pd_cairo)
    fmt_group.add_argument("--pil-png", action='store_const',
                           dest="draw_pd",
                           const=draw_pd_pil)
    parser.set_defaults(draw_pd=draw_pd_cairo)

    wh_size_group = parser.add_argument_group()
    wh_size_group.add_argument("-W", "--width", type=int,
                               help="Fixed width for resultant image (noop if scale)")
    wh_size_group.add_argument("-H", "--height", type=int,
                               help="Fixed height for resultant image (noop if scale)")
    parser.add_argument("-s", "--scale", type=float,
                        help="Fixed scale for resulting image grid")
    xyscale_size_group = parser.add_argument_group()
    xyscale_size_group.add_argument("-x", "--xscale", type=float,
                                    help="Fixed X scale for resulting image grid")
    xyscale_size_group.add_argument("-y", "--yscale", type=float,
                                    help="Fixed Y scale for resulting image grid")
    parser.set_defaults(width=100, height=100)

    parser.add_argument("-l", "--line_width", type=int,
                        default=1,
                        help="Width of the lines used to draw")
    parser.add_argument("-X", "--cross_width", type=int,
                        default=2,
                        help="Size of the gaps used to denote crossings")
    parser.add_argument("-r", "--cross_radius", type=int,
                        default=0,
                        help="Size of circle nodes used to denote crossings")

    subparsers = parser.add_subparsers()

    parser_pdstor = subparsers.add_parser(
        'pdstor', help='create images for diagrams in a pdstor')
    parser_pdstor.add_argument("pdstor", type=argparse.FileType('r'))
    parser_pdstor.add_argument("-o", "--skip-header", action='store_false',
                               dest="read_header")
    parser_pdstor.set_defaults(func=pdstor_command)

    parser_pdcode = subparsers.add_parser(
        'pdcode', help='create images from KnotTheory PDCodes')
    parser_pdcode.add_argument("pdcode", type=argparse.FileType('r'),
                               help="KnotTheory PDCode file to load")
    parser_pdcode.set_defaults(func=pdcode_command)

    parser_random = subparsers.add_parser(
        'random', help='create image from a random diagram')
    parser_random.add_argument("size", type=int)
    parser_random.set_defaults(func=random_command)

    args = parser.parse_args()
    print args
    args.func(args)

