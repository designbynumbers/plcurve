---------------------------------
README : Knot Shadow Data
---------------------------------

This directory contains supplementary data files for the paper "Knot
Probabilities in Random Diagrams" by Jason Cantarella, Harrison
Chapman, and Matt Mastin. This document describes the files and their
formats.

Due to the extremely tight filesize limits on arXiv ancillary data, we are
unable to distribute files of shadows or diagrams directly through arXiv.
They are currently available at

www.jasoncantarella.com/wordpress/papers

This is clearly not a stable solution, but we hope to enter the data in a
more permanent repository as the paper is published.

We also give csv files giving the actual counts of all observed knot
types for each number of crossings. Knot types are specified using a
version of Rolfsen's table augmented with composite knots and data on
the symmetry type of the knot. These tables appeared first in ``The
Shapes of Tight Composite Knots'' (Cantarella, LaPointe, and Rawdon)
and ``An Enhanced Prime Decomposition Theorem for Knots with Symmetry
Information'' (Mastin), and are included in this package as

knot-table.csv

In each, knot types are given in the format

Knot[<Crossing number>,<Index in Rolfsen's table>]<symmetry>

where <symmetry> is one of 'm' (all crossing signs reversed, or
'mirror image'), 'r' (orientation reversed), 'mr' (crossing signs and
orientation reversed), or blank. We took the versions of the knots in
the table of knots in the KnotTheory package to be the 'base' versions
of each knot. For composite knots, we give the names of the prime
factors with '#' in between.

The csv files of observations

knot-frequency-03.csv

through

knot-frequency-10.csv

are in the format

<knot type>, <number of diagrams>

with two special knot types:

Unknot[] and AllDiagrams[]

to indicate the number of unknots and the total number of diagrams at
the start of the file.




