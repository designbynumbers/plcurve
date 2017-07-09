#!/usr/bin/python
from libpl.pdcode import PlanarDiagram, UNSET_ORIENTATION
from libpl import pdcode

from collections import defaultdict

import json
from importlib import import_module

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
                          scale=None, xscale=None, yscale=None,
                          border=None, **kwargs):
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
        verts = [(v[0]-scale+border, v[1]-scale+border) for v in verts]
        maxx = max(v[0] for v in verts)
        maxy = max(v[1] for v in verts)
        width, height = border+maxx, border+maxy

    elif xscale is not None and yscale is not None:
        verts, arrows, xings = orth.ascii_data(emb, xscale, yscale)
        verts = [(v[0]-scalex+border, v[1]-scaley+border) for v in verts]
        maxx = max(v[0] for v in verts)
        maxy = max(v[1] for v in verts)
        width, height = border+maxx, border+maxy
    else:
        verts, arrows, xings = orth.ascii_data(emb, 1., 1.)

        maxx = max(v[0] for v in verts)
        maxy = max(v[1] for v in verts)

        scale = min(maxx, maxy)
        xscale = (width-2*border)/(maxx-1)
        yscale = (height-2*border)/(maxy-1)
        ssc = min(xscale, yscale)

        verts, arrows, xings = orth.ascii_data(emb, xscale, yscale)
        verts = [(v[0]-xscale+border, v[1]-yscale+border) for v in verts]

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
        self.cairo = import_module("cairo")
        self.fname = fname
        self.orient = True
        self.thickness = 1.

    def open(self, width, height):
        self.surf = self.cairo.SVGSurface(self.fname, width, height)
        self.cr = self.cairo.Context(self.surf)
        self.cr.set_line_width(self.thickness)

    def set_thickness(self, thickness):
        self.cr.set_line_width(thickness)
        self.thickness = thickness

    def set_color(self, r, g, b):
        self.cr.set_source_rgb(r, g, b)

    def draw_line(self, start, end):
        self.cr.move_to(*start)
        self.cr.line_to(*end)

        if self.orient:
            # Draw arrows
            cx = (end[0])
            cy = (end[1])

            if start[0] < end[0]:
                self.cr.move_to(cx-self.thickness*4, cy-self.thickness*3)
                self.cr.line_to(cx, cy)
                self.cr.line_to(cx-self.thickness*4, cy+self.thickness*3)
            elif start[0] > end[0]:
                self.cr.move_to(cx+self.thickness*4, cy+self.thickness*3)
                self.cr.line_to(cx, cy)
                self.cr.line_to(cx+self.thickness*4, cy-self.thickness*3)
            elif start[1] < end[1]:
                self.cr.move_to(cx+self.thickness*3, cy-self.thickness*4)
                self.cr.line_to(cx, cy)
                self.cr.line_to(cx-self.thickness*3, cy-self.thickness*4)
            elif start[1] > end[1]:
                self.cr.move_to(cx+self.thickness*3, cy+self.thickness*4)
                self.cr.line_to(cx, cy)
                self.cr.line_to(cx-self.thickness*3, cy+self.thickness*4)

        self.cr.stroke()


    def draw_circle(self, center, radius):
        self.cr.arc(center[0], center[1], radius, 0, 2*pi)
        self.cr.fill()

    def close(self):
        pass

class PILDrawingContext(DrawingContext):
    def __init__(self, fname, fmt="GIF"):
        self.PILImage = import_module("PIL.Image")
        self.PILImageDraw = import_module("PIL.ImageDraw")
        self.fname = fname
        self.color = (0,0,0)
        self.thickness = 1
        self.fmt = fmt

    def open(self, width, height):
        self.image = self.PILImage.new("RGB", (int(width), int(height)), (255,255,255))
        self.draw = self.PILImageDraw.Draw(self.image)

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
        with open(self.fname, "w") as f:
            self.image.save(f, self.fmt)

def _draw_pd(pd, dc, skip_start=0, skip_end=0,
             cross_width=4, cross_radius=0,
             line_width=1, **kwargs):
    orth = OrthogonalPlanarDiagram(pd)
    width, height, easy_arrows, broken_arrows = get_orthogonal_arrows(
        orth, **kwargs)

    dc.open(width, height)
    #print width,height

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

def draw_pd_pil(pd, fname, draw_pd, **kwargs):
    _draw_pd(pd, PILDrawingContext(fname, fmt=draw_pd["fmt"]),
             **kwargs)

from math import pi

def _args_to_draw_kwargs(args):
    return {"height": args.height,
            "width": args.width,
            "scale": args.scale,
            "xscale": args.xscale,
            "yscale": args.yscale,
            "border": args.border,
            "line_width": args.line_width,
            "cross_width": args.cross_width,
            "cross_radius": args.cross_radius,
            "draw_pd": args.draw_pd}

def filename_gen(uid, fmt, prefix=None):
    if prefix is not None:
        return "{}_{}.{}".format(prefix, uid, fmt)
    else:
        return "{}.{}".format(uid, fmt)

def pdstor_command(args):
    #print args.pdstor.name
    prefix = None
    if args.prefix:
        prefix = args.prefix
    elif args.auto_prefix:
        import os.path
        if args.pdstor.name[0] != "<":
            # Not a special file descriptor
            prefix = os.path.splitext(os.path.basename(args.pdstor.name))[0]

    gen_images = []
    for pd in PlanarDiagram.read_all(args.pdstor,
                                     read_header=args.read_header):
        if any(x.sign == UNSET_ORIENTATION for x in pd.crossings):
            # We don't know how to draw unset crossing signs yet...
            pd.randomly_assign_crossings()
        fname = filename_gen(pd.uid, args.draw_pd["ext"], prefix)
        if args.gallery:
            gen_images.append(fname)
        args.draw_pd["func"](pd,
                             fname,
                             **_args_to_draw_kwargs(args))

    MAX_GALLERY = 1000
    if args.gallery > 0:
        if len(gen_images) < MAX_GALLERY:
            with open("gallery.html", "w") as f:
                f.write("<html><body>\n")
                for i, src in enumerate(gen_images):
                    if i and not i%args.gallery:
                        f.write("<br>\n")
                    f.write("<img src='{}'>".format(src))
                f.write("</body></html>\n")
        else:
            gal = 0
            while gal*MAX_GALLERY < len(gen_images):
                with open("gallery_{}.html".format(gal), "w") as f:
                    f.write("<html><body>\n")
                    for i, src in enumerate(
                            gen_images[gal*MAX_GALLERY:gal*MAX_GALLERY+MAX_GALLERY+1]):
                        if i and not i%args.gallery:
                            f.write("<br>\n")
                        f.write("<img src='{}'>".format(src))
                    f.write("</body></html>\n")
                    gal += 1

def pdcode_command(args):
    i = 0
    pd = PlanarDiagram.read_knot_theory(args.pdcode)
    while pd is not None:
        args.draw_pd(pd,
                     filename_gen(i, args.draw_pd["ext"], args.prefix),
                     **_args_to_draw_kwargs(args))
        pd = PlanarDiagram.read_knot_theory(args.pdcode)
        i += 1

def random_command(args):
    pd = PlanarDiagram.random_diagram(args.size)
    args.draw_pd(pd,
                 filename_gen(pd.uid, args.draw_pd["ext"], args.prefix),
                 **_args_to_draw_kwargs(args))

format_opts = {
    'cairo-svg': {'func': draw_pd_cairo, "ext": "svg"},
    'pil-png': {"func": draw_pd_pil, "fmt": "PNG", "ext": "png"},
    'pil-gif': {"func": draw_pd_pil, "fmt": "GIF", "ext": "gif"}
}

if __name__ == "__main__":
    import sys, argparse

    class DictChoiceAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            super(DictChoiceAction, self).__init__(option_strings, dest, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, self.choices[values])

    parser = argparse.ArgumentParser(
        description='Generate SVG images from PlanarDiagrams')

    parser.add_argument("-f", "--format", type=str, choices=format_opts,
                        action=DictChoiceAction, dest="draw_pd")
    # fmt_group = parser.add_mutually_exclusive_group()
    # fmt_group.add_argument("--cairo-svg", action='store_const',
    #                        dest="draw_pd",
    #                        const={"func": draw_pd_cairo, "ext": "svg"})
    # fmt_group.add_argument("--pil-png", action='store_const',
    #                        dest="draw_pd",
    #                        const=format_opts['cairo-svg'])
    parser.set_defaults(draw_pd={"func": draw_pd_cairo, "ext": "svg"})

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

    parser.add_argument("-b", "--border", type=int,
                        help="Size of border around image",
                        default=5)

    parser.add_argument("-l", "--line_width", type=int,
                        default=1,
                        help="Width of the lines used to draw")
    parser.add_argument("-X", "--cross_width", type=int,
                        default=2,
                        help="Size of the gaps used to denote crossings")
    parser.add_argument("-r", "--cross_radius", type=int,
                        default=0,
                        help="Size of circle nodes used to denote crossings")


    parser.add_argument("-p", "--prefix", type=str,
                        help="Prefix for image filenames: [prefix]_[uid].[fmt]")

    subparsers = parser.add_subparsers()

    parser_pdstor = subparsers.add_parser(
        'pdstor', help='create images for diagrams in a pdstor')
    parser_pdstor.add_argument("pdstor", type=argparse.FileType('r'))
    parser_pdstor.add_argument("-o", "--skip-header", action='store_false',
                               dest="read_header",
                               help="Don't read pdstor header (e.g. for naked pd output)")
    parser_pdstor.add_argument("-N", "--no-auto-prefix", action="store_false",
                               dest="auto_prefix",
                               help="Don't automagically generate a prefix from filename")
    parser_pdstor.add_argument("-G", "--gallery", type=int,
                               help="Auto-generate HTML gallery of GALLERY columns")
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
    #print args
    args.func(args)

