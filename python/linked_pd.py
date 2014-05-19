import libplcurve.pd as pd

class LinkedDiagram(object):
    class LinkedComponent(object):

        class LinkedEdge(object):
            def __init__(self):
                self.head, self.headpos = None, None
                self.tail, self.tailpos = None, None
            def attach_head(self, cross, pos):
                self.head = cross
                self.headpos = pos
                cross[pos] = self
            def attach_tail(self, cross, pos):
                self.tail = cross
                self.tailpos = pos
                cross[pos] = self

        class LinkedCrossing(object):
            def __init__(self, edges=None, sign=0):
                if edges is None:
                    self._edges = [None, None, None, None]
                else:
                    self._edges = edges
                self.sign = sign
            def __getitem__(self, i):
                return self._edges[i%4]
            def __setitem__(self, i, val):
                self._edges[i%4] = val

        def __init__(self, dia, comp):
            pd = dia._pd
            self.edges = {}
            new_edge = None
            new_crossing = None
            new_inpos = None

            j_0 = comp[0]
            for j in comp:
                edge = pd.edges[j]
                new_edge = self.LinkedEdge()
                if new_crossing:
                    new_edge.attach_tail(new_crossing, (new_inpos+2)%4)
                self.edges[j] = new_edge
                i = edge.head

                cross = pd.crossings[i]
                new_crossing = None
                new_inpos = edge.headpos
                if edge.head not in dia.crossings:
                    new_crossing = self.LinkedCrossing(sign=cross.sign)
                    dia.crossings[i] = new_crossing
                else:
                    new_crossing = dia.crossings[i]
                new_edge.attach_head(new_crossing, new_inpos)
            self.edges[j_0].attach_tail(new_crossing, (new_inpos+2)%4)

    def __init__(self, pd):
        self._pd = pd
        self.refresh_components()

    def refresh_components(self):
        self.crossings = {}
        self.components = [LinkedDiagram.LinkedComponent(self, comp)
                           for comp in self._pd.components]

dia = pd.PlanarDiagram.read(open("../8.ex.pdcode"))
link = LinkedDiagram(dia)
