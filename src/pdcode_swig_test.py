#!/usr/bin/python

import pdcode

pd = pdcode.Shadow.from_twist_knot(5)
print pd.component_and_pos(1)


pdA = chr(3)+chr(1)+chr(4)+chr(0)+chr(1)+chr(5)+chr(2)+chr(4)+chr(5)+chr(3)+chr(0)+chr(2)

pdB = chr(3)+chr(1)+chr(4)+chr(0)+chr(1)+chr(5)+chr(2)+chr(4)+chr(0)+chr(2)+chr(5)+chr(3)

pdC = chr(4)+chr(0)+chr(3)+chr(1)+chr(1)+chr(5)+chr(2)+chr(4)+chr(0)+chr(2)+chr(5)+chr(3)

pdD = chr(0)+chr(3)+chr(1)+chr(4)+chr(1)+chr(5)+chr(2)+chr(4)+chr(0)+chr(2)+chr(5)+chr(3)

print pdcode.isomorphic(pdA,3,pdC,3)

print pdA
