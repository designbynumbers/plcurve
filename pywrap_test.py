#!/usr/bin/python

# Importing should work
import libplcurve.plcurve as pl
import pickle

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

# Free the plCurve. It should no longer take up memory
del plc

# Create a new random closed polygon
plc = pl.PlCurve.random_closed_polygon(r, 425)

# Classify the polygon
kt, nposs = plc.classify()
print kt, nposs

print kt.num_factors

kt_str = "#".join("{}_{}".format(f.cr, f.ind) for f in kt.factors)
print kt_str

## Cleanup. Memory usage should be back to 0.
del plc
del kt
del r

# Create a gsl random number generator
r = pl.RandomGenerator()

# Set a seed of 5000
r.set(5000)

# Generate a random planar PlCurve
planar = pl.PlCurve.random_closed_plane_polygon(r, 10)

# Generate another planar PlCurve
planar2 = pl.PlCurve.random_closed_plane_polygon(r, 10)

# Ensure the two are distinct (probability 0 they are equal...)
assert(planar.components[0].vertices != planar2.components[0].vertices)

# Create a new RNG from the state of the old
r2 = pl.RandomGenerator.from_state(r.get_state())

# The new and old RNG should now generate the same random polygon
poly1 = pl.PlCurve.random_closed_polygon(r, 10)
poly2 = pl.PlCurve.random_closed_polygon(r2, 10)
assert(poly1.components[0].vertices == poly2.components[0].vertices)

### Demonstrate file i/o for plCurves
# Load a plcurve from data
f = open("data/knotvects/9_40.vect", "r")
plc = pl.PlCurve.from_file(f)
f.close()

#print plc.components[0].vertices
# Classify the knot to check that it is as it purports
kt, nposs = plc.classify()
kt_str = "#".join("{}_{}".format(f.cr, f.ind) for f in kt.factors)
print kt_str

# should work, but doesn't: TODO: Fix!!
#print "{}_{}".format(str(kt.factors.cr), str(kt.factors.ind))

# Write a random plcurve to a file
plc = pl.PlCurve.random_closed_polygon(r, 40)
f = open("____definitely_fake_plcurve.vect", "w")
plc.write(f)
f.close()
# file should now be written
