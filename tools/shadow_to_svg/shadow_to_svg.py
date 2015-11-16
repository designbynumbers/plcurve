import svgwrite
from spherogram.links.orthogonal import OrthogonalLinkDiagram

def shadow_to_svg(pd, fname="shadow.svg"):
    orth_dia = OrthogonalLinkDiagram(pd.as_spherogram())
    verts, arrows, xings = orth_dia.plink_data()
    arrows = dict(arrows)

    components = []
    comp_start = None
    arrw_start = comp_start
    component = []
    while arrows:
        if arrw_start is None:
            arrow = arrows.popitem()
            comp_start = arrow[0]
        else:
            arrow = arrw_start, arrows[arrw_start]
            del arrows[arrw_start]

        arrw_start = arrow[1]
        component.append(arrow[0])

        if arrw_start == comp_start:
            components.append(component)
            component = []
            arrw_start = None
            comp_start = None

    assert(component == [])

    SCALE=4

    xs = [x*SCALE/10 for x,y in verts]
    ys = [y*SCALE/10 for x,y in verts]

    max_x = max(xs)+SCALE
    max_y = max(ys)+SCALE

    dwg = svgwrite.Drawing(size=(max_x,max_y))
    for component in components:
        comp_verts = [verts[i] for i in component]
        comp_verts = [(x*SCALE/10,y*SCALE/10) for x,y in comp_verts]
        poly = dwg.add(dwg.polygon(comp_verts,
                                   fill='none',
                                   stroke='black',
                                   stroke_width=.4*SCALE))

    dwg.saveas(fname)
