plcurve: Polygonal Links
------------------------

.. currentmodule:: libplcurve.plcurve

The PlCurve class
'''''''''''''''''

.. autoclass:: PlCurve

   In addition to the methods below, PlCurves have the following
   functionality:

   + A PlCurve object can act like a list of its components

      + Access of component by index: ::

          plc = pl.PlCurve.from_file("example_pl.vert")
          first_component = plc[0]
          also_first_component = plc.components[0]

      + Deletion of component by index with ``del``: ::

          del plc[0]

      + The length of a PlCurve is its number of vertices: ::

          len(plc)

      + Iteration through components works as expected: ::

          for component in plc:
              my_fancy_component_op(component)
              print component

   + A PlCurve has a string representation which briefly describes it: ::

       >>> print plc
       PlCurve with 3 components

   **PlCurve constructor**

   **Static methods which create new PlCurves**

   .. automethod:: from_file

   **Static methods which create new random PlCurves**

   These functions use the gsl random number generator facilities.
   This adds some (minor) responsibilities for user programs; namely,
   you must create a random number generator explicitly before calling
   these functions.

   Typically, this will involve creating a :py:class:`RandomGenerator`
   instance and (optionally) setting a seed with
   :py:meth:`RandomGenerator.set`: ::

       import libplcurve.plcurve as pl
       random_gen = pl.RandomGenerator()
       random_gen.set(12345)

       random_plc = pl.PlCurve.random_closed_polygon(r, 50)

   .. automethod:: random_closed_polygon

   .. automethod:: random_open_polygon

   .. automethod:: random_closed_plane_polygon

   .. automethod:: random_open_plane_polygon

   .. automethod:: random_equilateral_closed_polygon

   .. automethod:: random_equilateral_open_polygon

   .. automethod:: loop_closure

   **Data operation methods**

   .. automethod:: add_component

   .. automethod:: append

   .. automethod:: drop_component

   .. automethod:: resize_colorbuf

   .. automethod:: set_fixed

   .. automethod:: constrain_to_line

   .. automethod:: constrain_to_plane

   .. automethod:: unconstrain

   .. automethod:: remove_constraint

   .. automethod:: remove_all_constraints

   .. automethod:: fix_wrap

   .. automethod:: fix_constraints

   .. automethod:: double_vertices

   **Spline methods**

   .. automethod:: convert_to_spline

   **Geometric information method & attributes**

   .. automethod:: num_edges

   .. automethod:: num_vertices

   .. automethod:: vertex_num

   .. automethod:: cp_num

   .. automethod:: vt_num

   .. automethod:: turning_angle

   .. automethod:: MR_curvature

   .. automethod:: total_curvature

   .. automethod:: total_torsion

   .. automethod:: mean_tangent

   .. automethod:: subarc_length

   .. automethod:: s

   .. automethod:: check_constraint

   .. automethod:: mean_squared_chordlengths

   .. autoattribute:: pointset_diameter

      Calculate the diameter of the PlCurve, thinking of vertices as set of
      points in \\(\\mathbb{R}^3\\).

   .. autoattribute:: center_of_mass

   .. autoattribute:: gyradius

   **Geometric operation methods**

   .. automethod:: scale

   .. automethod:: fold_move

   .. automethod:: rotate

   .. automethod:: random_rotate

   .. automethod:: translate

   .. automethod:: perturb

   .. automethod:: project

   .. automethod:: delete_arc

   **Symmetry methods**

   (Not yet implemented)

   **Topology methods**

   .. automethod:: classify

   .. automethod:: homfly

   .. automethod:: as_pd

The Strand class
''''''''''''''''

.. autoclass:: Strand
   :members:

The RandomGenerator (gsl_rng interface) class
'''''''''''''''''''''''''''''''''''''''''''''

.. autoclass:: RandomGenerator
   :members:

Auxiliary classes
'''''''''''''''''

.. autoclass:: Symmetry
   :members:

.. autoclass:: SymmetryGroup
   :members:

.. autoclass:: Symmetry
   :members:

.. autoclass:: Constraint
   :members:

.. autoclass:: VertexQuantifier
   :members:

.. autoclass:: SplineStrand
   :members:

.. autoclass:: Spline
   :members:

.. autoclass:: KnotType
   :members:
