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

The original GDSII format limits the number of vertices in a polygon to 199.
Most modern software disregards this limit and allows an arbitrary number of points per polygon.
Gdspy follows the modern version of GDSII, but this is an important issue to keep in mind if the generated file is to be used in older systems.

The units used to represent shapes in the GDSII format are defined by the user.
The default unit in gdspy is 1 µm (10⁻⁶ m), but that can be easily changed by the user.


***********
First GDSII
***********

Let's create our first GDSII file:

.. code-block:: python

   import gdspy

   # Create the geometry: a single rectangle.
   rect = gdspy.Rectangle((0, 0), (2, 1))
   cell = gdspy.Cell('FIRST')
   cell.add(rect)

   # Save all created cells in file 'first.gds'.
   gdspy.write_gds('first.gds')

   # Optionally, display all cells using the internal viewer.
   gdspy.LayoutViewer()


After importing the gdspy module, we create a :class:`gdspy.Rectangle` with opposing corners at positions (0, 0) and (2, 1).

Then a :class:`gdspy.Cell` is created and the rectangle is added to the cell.
All shapes in the GDSII format exist inside cells.
A cell can be imagined as a piece of paper where the layout will be defined.
Later, the cells can be used to create a hierarchy of geometries, ass we'll see in :ref:`References`.

Finally, the whole structure is saved in a file called "first.gds" in the current directory.
By default, all created cells are included in this operation.

The GDSII file can be opened in a number of viewers and editors, such as `KLayout <https://klayout.de/>`_.
Alternatively, gdspy includes a simple viewer that can also be used: :class:`gdspy.LayoutViewer`.


********
Polygons
********

General polygons can be defined by an ordered list of vertices.
The orientation of the vertices (clockwise/counter-clockwise) is not important.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Polygons
   :end-before: draw

.. image:: _static/polygons.*
   :align: center


Holes
=====

As mentioned in :ref:`Getting Started`, holes have to be connected to the outer boundary of the polygon, as in the following example:

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Holes
   :end-before: draw

.. image:: _static/holes.*
   :align: center


Circles
=======

The :class:`gdspy.Round` class creates circles, ellipses, doughnuts, arcs and slices.
In all cases, the arguments `tolerance` or `number_of_points` will control the number of vertices used to approximate the curved shapes.

If the number of vertices in the polygon is larger than `max_points` (199 by default), it will be fractured in many smaller polygons with at most `max_points` vertices each.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Circles
   :end-before: draw

.. image:: _static/circles.*
   :align: center


Transformations
===============

All polygons can be transformed trough :meth:`gdspy.PolygonSet.translate`, :meth:`gdspy.PolygonSet.rotate`, :meth:`gdspy.PolygonSet.scale`, and :meth:`gdspy.PolygonSet.mirror`.
The transformations are applied in-place, i.e., no polygons are created.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Transformations
   :end-before: draw

.. image:: _static/transformations.*
   :align: center


Layer and Datatype
==================

All shapes in the GDSII format are tagged with 2 properties: layer and datatype (or texttype in the case of :class:`gdspy.Label`).
They are always 0 by default, but can be any integer in the range from 0 to 255.

These properties have no predefined meaning.
It is up to the system using the GDSII file to chose with to do with those tags.
For example, in the CMOS fabrication process, each layer could represent a different lithography level.

In the example below, a single file stores different fabrication masks in separate layer and datatype configurations.


.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Layer and Datatype
   :end-before: draw

.. image:: _static/layer_and_datatype.*
   :align: center


**********
References
**********

References give the GDSII format its hierarchical features.
They work by reusing a cell content in another cell (without actually copying the whole geometry).
As a simplistic example, imagine the we are designing a simple electronic circuit that uses hundreds of transistors, but they all have the same shape.
We can draw the transistor just once and reference it throughout the circuit, rotating or mirroring each instance as necessary.

Besides creating single references with :class:`gdspy.CellReference`, it is possible to create full 2D arrays with a single entity using :class:`gdspy.CellArray`.
Both are exemplified below.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: References
   :end-before: draw

.. image:: _static/references.*
   :align: center


*****
Paths
*****

Besides polygons, the GDSII format defines paths, witch are `polygonal chains <https://en.wikipedia.org/wiki/Polygonal_chain>`_ with associated width and end caps.
The width is a single number, constant throughout the path, and the end caps can be flush, round, or extended by a custom distance.

There is no specification for the joins between adjacent segments, so it is up to the system using the GDSII file to specify those.
Usually the joins are straight extensions of the path boundaries up to some beveling limit.
Gdspy also uses this specification for the joins.

It is possible to circumvent all of the above limitations within gdspy by storing paths as polygons in the GDSII file.
The disadvantage of this solution is that other software will not be able to edit the geometry as paths, since that information is lost.

The construction of paths (either GDSII paths or polygonal paths) in gdspy is quite rich.
There are 3 classes that can be used depending on the requirements of the desired path.


Polygonal-Only Paths
====================

The class :class:`gdspy.Path` is designed to allow the creation of path-like polygons in a piece-wise manner.
It is the most computationally efficient class between the three because it *does not* calculate joins.
That means the user is responsible for designing the joins.
The paths can end up with discontinuities if care is not taken when creating them.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Polygonal-Only Paths
   :end-before: draw

.. image:: _static/polygonal-only_paths.*
   :align: center

Just as with :ref:`Circles`, all curved geometry is approximated by line segments.
The number of segments is similarly controlled by a `tolerance` or a `number_of_points` argument.
Curves also include fracturing to limit the number of points in each polygon.

More complex paths can be constructed with the methods :meth:`gdspy.Path.bezier`, :meth:`gdspy.Path.smooth`, and :meth:`gdspy.Path.parametric`.
The example below demonstrates a couple of possibilities.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Polygonal-Only Paths 1
   :end-before: draw

.. image:: _static/polygonal-only_paths_1.*
   :align: center

The width of the path does not have to be constant.
Each path component can linearly taper the width of the path by using the `final_width` argument.
In the case of a parametric curve, more complex width changes can be created by setting `final_width` to a function.

Finally, parallel paths can be created simultaneously with the help of arguments `number_of_paths`, `distance`, and `final_distance`.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Polygonal-Only Paths 2
   :end-before: draw

.. image:: _static/polygonal-only_paths_2.*
   :align: center


Simple Paths
============

Although very efficient, :class:`gdspy.Path` is limited in the type of path it can provide.
For example, if we simply want a path going through a sequence of points, we need a class that can correctly compute the joins between segments.
That's one of the advantages of class :class:`gdspy.SimplePath`.
Other path construction methods are similar to those in :class:`gdspy.Path`.

A few features of :class:`gdspy.SimplePath` are:

- paths can be stored as proper GDSII paths;
- end caps and joins can be specified by the user;
- each parallel path can have a different width;
- spacing between parallel paths is arbitrary; the user specifies the offset of each path individually.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Simple Paths
   :end-before: draw

.. image:: _static/simple_paths.*
   :align: center

The following example shows other features, such as width tapering, arbitrary offsets, and custom joins and end caps.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Simple Paths 1
   :end-before: draw

.. image:: _static/simple_paths_1.*
   :align: center


The corner type 'circular bend' (together with the `bend_radius` argument) can be used to automatically curve the path.
This feature is used in :ref:`Example: Integrated Photonics`.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Simple Paths 2
   :end-before: draw

.. image:: _static/simple_paths_2.*
   :align: center


Lazy Paths
==========

In some situations, :class:`gdspy.SimplePath` is unable to properly calculate all the joins.
This often happens when the width or offset of the path is relatively large with respect to the length of the segments being joined.
Curves that join other curves or segments at sharp angles are an example of such situation.

The class :class:`gdspy.LazyPath` can be used in such scenarios where robustness is more important than efficiency due to sharp corners or large offsets in the paths.
The drawbacks of using :class:`gdspy.LazyPath` are the loss in computation efficiency (compared to the other 2 classes) and the impossibility of specifying corner shapes.
The advantages are, as mentioned earlier, more robustness when generating the final geometry, and freedom to use custom functions to parameterize the widths or offsets of the paths in any construction method.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Lazy Paths
   :end-before: draw

.. image:: _static/lazy_paths.*
   :align: center

Note that, analogously to :class:`gdspy.SimplePath`, :class:`gdspy.LazyPath` can be stored as a GDSII path as long as its width is kept constant.


****
Text
****

In the context of a GDSII file, text is supported in the form of labels, which are ASCII annotations placed somewhere in the geometry of a given cell.
Similar to polygons, labels are tagged with layer and texttype values (texttype is the label equivalent of the polygon datatype).
They are supported by the class :class:`gdspy.Label`.

Additionally, gdspy offers the possibility of creating text as polygons to be included with the geometry.
The class :class:`gdspy.Text` creates polygonal text that can be used in the same way as any other polygons in gdspy.
The font used to render the characters contains only horizontal and vertical edges, which is important for some laser writing systems.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: # Text
   :end-before: draw

.. image:: _static/text.*
   :align: center


*******************
Geometry Operations
*******************

Gdspy offers a number of functions and methods to modify existing geometry.
The most useful operations include :func:`gdspy.boolean`, :func:`gdspy.slice`, :func:`gdspy.offset`, and :meth:`gdspy.PolygonSet.fillet`.


Boolean Operations
==================

Boolean operations (:func:`gdspy.boolean`) can be performed on polygons, paths and whole cells.
Four operations are defined: union ('or'), intersection ('and'), subtraction ('not'), and symmetric subtraction ('xor').

They can be computationally expensive, so it is usually advisable to avoid using boolean operations whenever possible.
If they are necessary, keeping the number of vertices is all polygons as low as possible also helps.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Boolean Operations
   :end-before: draw

.. image:: _static/boolean_operations.*
   :align: center


Slice Operation
===============

As the name indicates, a slice operation subdivides a set of polygons along horizontal or vertical cut lines.

In a few cases, a boolean operation can be substituted by one or more slice operations.
Because :func:`gdspy.slice` is ususally much simpler than :func:`gdspy.boolean`, it is a good idea to use the former if possible.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Slice Operation
   :end-before: draw

.. image:: _static/slice_operation.*
   :align: center


Offset Operation
================

The function :func:`gdspy.offset` expands or contracts polygons by a fixed amount.
It can operate on individual polygons or sets of them, in which case it may make sense to use the argument `join_first` to operate on the whole geometry as if a boolean 'or' was executed beforehand.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Offset Operation
   :end-before: draw

.. image:: _static/offset_operation.*
   :align: center


Fillet Operation
================

The method :meth:`gdspy.PolygonSet.fillet` can be used to round polygon corners.
It doesn't have a `join_first` argument as :func:`gdspy.offset`, so if it will be used on a polygon, that polygon should probably not be fractured.

.. literalinclude:: makeimages.py
   :language: python
   :dedent: 4
   :start-after: Fillet Operation
   :end-before: draw

.. image:: _static/fillet_operation.*
   :align: center


*************
GDSII Library
*************

TODO: Units, precision, `current_library`


Saving a GDSII File
===================


Loading a GDSII File
====================


*****************************
Example: Integrated Photonics
*****************************

This example demonstrates the use of gdspy primitives to create more complex structures.

These structures are commonly used in the field of integrated photonics.

:download:`photonics.py <_static/photonics.py>`

.. literalinclude:: _static/photonics.py
   :language: python
   :linenos:
