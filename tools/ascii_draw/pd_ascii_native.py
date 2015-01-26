#!/usr/bin/python

import pd_ascii
from libpl.pdcode import *

K = PlanarDiagram.torus_knot(2,3)
#K = PlanarDiagram.unknot(1)
O = pd_ascii.OrthogonalPlanarDiagram(K)
pd_ascii.draw(O)
