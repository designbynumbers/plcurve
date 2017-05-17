The Diagram module: PlanarDiagrams
==================================

.. automodule:: libpl.pdcode.diagram

Constants and other top-level features
--------------------------------------

The Diagram module provides the following constant values:

**Signs and orientation constants**

.. data:: POS_ORIENTATION

   Positive orientation sign

.. data:: NEG_ORIENTATION

   Negative orientation sign

.. data:: UNSET_ORIENTATION

   Unset orientation sign

**Equivalence type constants**

.. data:: ISOTOPY

   Use diagram isomorphism; check signs

.. data:: UO_ISOMORPHISM

   Use map iosmorphism; ignore signs

Additionally, the diagram module includes some helper functions for toggling
debug verbosity in the wrapped C library:

.. autofunction:: pd_debug_off
.. autofunction:: pd_debug_on

The PlanarDiagram Class
-----------------------

The primary definition of the module `diagram` is the `PlanarDiagram` object:

.. autoclass:: PlanarDiagram
   :members: ncomps, nedges, ncross, nfaces, hash, uid, edges, crossings,
             components, faces, thin, hashed, equiv_type

Creating PlanarDiagrams
'''''''''''''''''''''''

Typically, one does not construct such an object by hand. There are several
different ways to either create or load PlanarDiagram objects.

The following are all class methods:

.. automethod:: PlanarDiagram.random_diagram

.. automethod:: PlanarDiagram.from_dt_code
.. automethod:: PlanarDiagram.from_editor
.. automethod:: PlanarDiagram.from_plink
.. automethod:: PlanarDiagram.from_pdcode

.. automethod:: PlanarDiagram.db_knot
.. automethod:: PlanarDiagram.db_link

.. automethod:: PlanarDiagram.unknot
.. automethod:: PlanarDiagram.unknot_wye
.. automethod:: PlanarDiagram.torus_knot
.. automethod:: PlanarDiagram.twist_knot
.. automethod:: PlanarDiagram.simple_chain

The following methods produce new diagrams from other diagrams:

.. automethod:: PlanarDiagram.copy
.. automethod:: PlanarDiagram.edit_copy

Diagram I/O
'''''''''''

There are a few methods for reading and writing specific diagrams from and to
files. The following "read" methods are all class methods:

.. automethod:: PlanarDiagram.read
.. automethod:: PlanarDiagram.read_all
.. automethod:: PlanarDiagram.read_knot_theory

The following "write" methods are also available to save diagrams to disk.

.. automethod:: PlanarDiagram.write
.. automethod:: PlanarDiagram.write_knot_theory

The following method checks if, and where, a diagram is in a pdstor file.

.. automethod:: PlanarDiagram.pdstor_index

Diagram properties and invariants
'''''''''''''''''''''''''''''''''

It is possible to find properties and other invariants of diagrams and shadows.

.. automethod:: PlanarDiagram.homfly
.. automethod:: PlanarDiagram.ccode

.. automethod:: PlanarDiagram.linking_number
.. automethod:: PlanarDiagram.unsigned_linking_number
.. automethod:: PlanarDiagram.interlaced_crossings
.. automethod:: PlanarDiagram.interlaced_crossings_unsigned

.. automethod:: PlanarDiagram.euler_characteristic

Equivalence Checking and Isomorphism generation
'''''''''''''''''''''''''''''''''''''''''''''''

In addition to the ``==`` and ``!=`` comparisons (which depend on the value of
`equiv_type`), diagrams support the following operations for comparison and
equivalence checking:

.. automethod:: PlanarDiagram.isotopic
.. automethod:: PlanarDiagram.map_isomorphic
.. automethod:: PlanarDiagram.isomorphic
.. automethod:: PlanarDiagram.identical

Furthermore, it is possible to inspect isomorphism groups using,

.. automethod:: PlanarDiagram.build_isotopies
.. automethod:: PlanarDiagram.build_autoisotopies
.. automethod:: PlanarDiagram.build_map_isomorphisms
.. automethod:: PlanarDiagram.build_isomorphisms
.. automethod:: PlanarDiagram.build_automorphisms

Diagram Manipulation
''''''''''''''''''''

There are several different methods for the manipulation of PlanarDiagram objects.

.. automethod:: PlanarDiagram.set_all_crossing_signs
.. automethod:: PlanarDiagram.unset_crossing_signs
.. automethod:: PlanarDiagram.toggle_crossing
.. automethod:: PlanarDiagram.reorient_component
.. automethod:: PlanarDiagram.randomly_assign_crossings

.. automethod:: PlanarDiagram.connect_sum

The following methods change diagrams at a lower level, and require some care.

.. automethod:: PlanarDiagram.resize
.. automethod:: PlanarDiagram.add_crossing

Simplification and Reidemeister Moves
'''''''''''''''''''''''''''''''''''''

One goal of the PlanarDiagram library is to implement routines for
simplification.

.. automethod:: PlanarDiagram.neighbors
.. automethod:: PlanarDiagram.simplify

This works (theoretically) by applying Reidemeister moves to diagrams which do
not change the link type of a diagram. The following methods check for
Reidemeister local sites:

.. automethod:: PlanarDiagram.get_R1_loops
.. automethod:: PlanarDiagram.get_R2_bigons
.. automethod:: PlanarDiagram.get_R3_triangles

The following Reidemeister-type moves are implemented:

.. automethod:: PlanarDiagram.R1_loop_deletion
.. automethod:: PlanarDiagram.R2_bigon_elimination
.. automethod:: PlanarDiagram.R2_bigon_elimination_vertices

The following Reidemeister methods are non-working, or otherwise non-tested

.. automethod:: PlanarDiagram.R1_loop_addition
.. automethod:: PlanarDiagram.R2_bigon_addition
.. automethod:: PlanarDiagram.R3_triangle_flip
.. automethod:: PlanarDiagram.R3_nice_swap
.. automethod:: PlanarDiagram.R3_strand_swap

Regeneration and Sanity Checking
''''''''''''''''''''''''''''''''

Sometimes---particularly after making low-level modifications to diagram
objects---we have to regenerate higher-level objects like :obj:`Component` and
:obj:`Face` and check that planar diagram objects are valid.

*Note*: Methods of the form ``regenerate_componentType`` do not regenerate
 Python objects, and should be used with care; it is advised to instead use the
 `~PlanarDiagram.regenerate` method instead.

.. automethod:: PlanarDiagram.regenerate
.. automethod:: PlanarDiagram.regenerate_hash

.. automethod:: PlanarDiagram.is_ok
.. automethod:: PlanarDiagram.crossings_ok
.. automethod:: PlanarDiagram.edges_ok
.. automethod:: PlanarDiagram.faces_ok
.. automethod:: PlanarDiagram.components_ok

.. automethod:: PlanarDiagram.assert_nonnull

SnapPy and Spherogram
'''''''''''''''''''''

We support some conversions to :py:mod:`snappy` and :py:mod:`spherogram` data
types.

.. automethod:: PlanarDiagram.as_spherogram
.. automethod:: PlanarDiagram.snappy_manifold
