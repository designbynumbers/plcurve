import svgwrite
from spherogram.links.orthogonal import OrthogonalLinkDiagram

def shadow_to_svg(pd):
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

    dwg = svgwrite.Drawing()
    for component in components:
        comp_verts = [verts[i] for i in component]
        comp_verts = [(x/10,y/10) for x,y in comp_verts]
        poly = dwg.add(dwg.polygon(comp_verts))

    dwg.saveas("shadow.svg")
