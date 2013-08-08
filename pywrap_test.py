#!/usr/bin/python

# Importing should work
import libplcurve.plcurve as pl

# Rudimentary way to create a RNG
r = pl.make_gsl_rng()

# Set up a randomly generated PlCurve
plc = pl.PlCurve.random_closed_polygon(r, 6)

# Should be a list of component objects
print plc.components

# Print out all of the vertices in the components
for comp in plc.components:
    print comp.vertices

# Print out a user-friendly string for the PlCurve
print(plc)
print(str(plc))
print(plc.__str__()) # This is the correct result
