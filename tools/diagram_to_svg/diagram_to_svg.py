#!/usr/bin/python
from libpl.pdcode import PlanarDiagram
from libpl import pdcode

from collections import defaultdict

import cairo
import json

from orthogonal import OrthogonalPlanarDiagram

with open("config.json") as f:
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

def draw_orthogonal_projection(orth, fname,
                               width=100., height=100.,
                               skip_start=0, skip_end=0,
                               cross_width=4):
    """Draw an OrthogonalPlanarDiagram"""

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

    surf = cairo.SVGSurface(fname, width, height)
    cr = cairo.Context(surf)
    for edge, arrowtail, arrowhead, comp in easy_arrows:
        # skip an easy arrow if its edge is skippable
        if skip_start <= edge <= skip_end:
            continue
        elif (skip_end < skip_start and
              (edge <= skip_end or edge >= skip_start)):
            continue

        cr.set_source_rgb(*colors[comp])
        cr.move_to(*arrowtail)
        cr.line_to(*arrowhead)
        cr.stroke()

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
                cr.set_source_rgb(*colors[comp])
                cr.move_to(arrowtail[0], arrowtail[1])
                cr.line_to(crossloc[0] - sgn*cross_width, crossloc[1])
                cr.move_to(crossloc[0] + sgn*cross_width, crossloc[1])
                cr.line_to(arrowhead[0], arrowhead[1])
                cr.stroke()
            else:
                sgn = 1 if arrowtail[1] < arrowhead[1] else -1
                #print "vtunder", arrowhead, arrowtail, crossloc
                cr.set_source_rgb(*colors[comp])
                cr.move_to(arrowtail[0], arrowtail[1])
                cr.line_to(crossloc[0], crossloc[1] - sgn*cross_width)
                cr.move_to(crossloc[0], crossloc[1] + sgn*cross_width)
                cr.line_to(arrowhead[0], arrowhead[1])
                cr.stroke()
        else:
            cr.set_source_rgb(*colors[comp])
            cr.move_to(*arrowtail)
            cr.line_to(*arrowhead)
            cr.stroke()
            cr.arc(crossloc[0], crossloc[1], cross_width/2, 0, 2*pi)
            cr.fill()

    surf.write_to_png("knot.png")

from math import pi

def pdstor_command(args):
    for pd in PlanarDiagram.read_all(args.pdstor, read_header=True):
        pd.randomly_assign_crossings()
        orth = OrthogonalPlanarDiagram(pd)
        draw_orthogonal_projection(orth, "diagram_{}.svg".format(pd.uid),
                                   height=400, width=400,cross_width=2)

def pdcode_command(args):
    for i in range(10):
        pd = PlanarDiagram.random_diagram(3, 1, 0)
        orth = OrthogonalPlanarDiagram(pd)
        for edge in pd.edges:
            assert edge.parent is not None

if __name__ == "__main__":
    import sys, argparse

    parser = argparse.ArgumentParser(
        description='Generate SVG images from PlanarDiagrams')

    subparsers = parser.add_subparsers()

    parser_pdstor = subparsers.add_parser(
        'pdstor', help='create images for diagrams in a pdstor')
    parser_pdstor.add_argument("pdstor", type=argparse.FileType('r'))
    parser_pdstor.set_defaults(func=pdstor_command)

    parser_pdcode = subparsers.add_parser(
        'pdcode', help='create image from a KnotTheory PDCode')
    parser_pdcode.add_argument("pdcode", type=str,
                               help="KnotTheory PDCode to load")
    parser_pdcode.set_defaults(func=pdcode_command)

    args = parser.parse_args()
    args.func(args)

