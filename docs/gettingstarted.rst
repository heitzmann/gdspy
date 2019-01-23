###############
Getting Started
###############

GDSII files contain a hierarchical representation of any polygonal geometry.
They are mainly used in the microelectronics industry for the design of mask layouts, but are also employed in other areas.

Because it is a hierarchical format, repeated structures, such as identical transistors, can be defined once and referenced multiple times in the layout, reducing the file size.

There is one important limitation in the GDSII format: it only supports `weakly simple polygons <https://en.wikipedia.org/wiki/Simple_polygon>`_, that is, polygons whose segments are allowed to intersect, but not cross.

In particular, curves and shapes with holes are *not* directly supported.
Holes can be defined, nonetheless, by connecting their boundary to the boundary of the enclosing shape.
In the case of curves, they must be approximated by a polygon.
The number of points in the polygonal approximation can be increased to better approximate the original curve up to some acceptable error.

The units used to represent shapes in the GDSII format are defined by the user.
The default unit in gdspy is 1 µm (10⁻⁶ m), but that can be easily changed by the user.


***********
First GDSII
***********

Let's create our first GDSII file:

.. code-block:: python
   :linenos:

   import gdspy

   # Create the geometry: a single rectangle.
   rect = gdspy.Rectangle((0, 0), (2, 1))
   cell = gdspy.Cell('MAIN')
   cell.add(rect)

   # Save all created cells in file 'first.gds'.
   gdspy.write_gds('first.gds')

   # Optionally, display all cells using the internal viewer.
   gdspy.LayoutViewer()


After importing the gdspy module, we create a rectangle with opposing corners at positions (0, 0) and (2, 1).

Then a cell is created and the rectangle is added to the cell
All shapes in the GDSII format exist inside cells.
A cell can be imagined as a piece of paper where the layout will be defined.
Later, the cells can be used to create a hierarchy of geometries, ass we'll see in :ref:`sec-references`.

Finally, the whole structure is saved in a file called "first.gds" in the current directory.
By default, all created cells are included in this operation.

The GDSII file can be opened in a number of viewers and editors, such as `KLayout <https://klayout.de/>`_.
Alternatively, gdspy includes a simple viewer that can also be used: :class:`gdspy.LayoutViewer`.


********
Polygons
********

TODO


Circles
=======

TODO


Transformations
===============

TODO


.. _sec-references:

**********
References
**********

TODO


*****
Paths
*****

TODO


*******************
Geometry Operations
*******************

TODO


******************
GDSII Manipulation
******************

TODO: Units, precision, loading/saving...


*****************************
Example: Integrated Photonics
*****************************

TODO: Use paths


This example demonstrates the use of gdspy primitives to create more complex structures.

These structures are commonly used in the field of integrated photonics.

:download:`photonics.py <../examples/photonics.py>`

.. literalinclude:: ../examples/photonics.py
   :language: python
   :linenos:
