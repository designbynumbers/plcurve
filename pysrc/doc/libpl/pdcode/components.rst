The Components module
=====================

.. automodule:: libpl.pdcode.components

Crossings: The primary building blocks
--------------------------------------

A simple example: Creating a trefoil diagram
''''''''''''''''''''''''''''''''''''''''''''

The following code excerpt demonstrates how you might use Crossing objects to
manually construct a `~.diagram.PlanarDiagram` object in Python:

.. code-block:: python

   from libpl.pdcode import PlanarDiagram
   print "Hello"

The Crossing class
''''''''''''''''''

.. autoclass:: Crossing
   :members:

Edges
-----

.. autoclass:: Edge
   :members:

Link components
---------------

.. autoclass:: Component
   :members:

Faces
-----

.. autoclass:: Face
   :members:
