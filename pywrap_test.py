#!/usr/bin/python

# Importing should work
import libplcurve.plcurve as pl

# Create a gsl random number generator
r = pl.RandomGenerator()

# Set up a randomly generated PlCurve
plc = pl.PlCurve.random_open_polygon(r, 6)

# Should be a list of component objects
print plc.components

# Print out all of the vertices in the components
for comp in plc.components:
    print comp.vertices
    print(comp.__str__())

# Print out a user-friendly string for the PlCurve
#  -builtin seems to screw this up. All three prints should be the same as last
print(plc)
print(str(plc))
print(plc.__str__()) # This is the correct result

# Add a component to the PlCurve. This tests a difficult `in` typemap.
plc.add_component(1, [(0,0,0), (1,1,1), (2,0,2)], True, None)
plc.add_component(2, [(0,0,0), (3,3,3), (3,0,0)], True, [(0,0,1,1)])

# Print the PlCurve again. It should now have 2 components.
print(plc.__str__())

del plc
plc = pl.PlCurve.random_closed_polygon(r, 425)

kt, nposs = plc.classify()
print kt, nposs

print kt.num_factors

print kt.factors.cr, kt.factors.ind, kt.factors.sym

## Cleanup
del plc
del kt
del r
