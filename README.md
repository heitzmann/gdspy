# GDSPY README

[![Boost Software License - Version 1.0](https://img.shields.io/github/license/heitzmann/gdspy.svg)](http://www.boost.org/LICENSE_1_0.txt)
[![Documentation Status](https://readthedocs.org/projects/gdspy/badge/?version=stable)](https://gdspy.readthedocs.io/en/stable/?badge=stable)
[![tests](https://github.com/heitzmann/gdspy/actions/workflows/testing.yml/badge.svg)](https://github.com/heitzmann/gdspy/actions/workflows/testing.yml)
[![Appveyor Status](https://ci.appveyor.com/api/projects/status/pr49a6bhxvbqwocy?svg=true)](https://ci.appveyor.com/project/heitzmann/gdspy)
[![Downloads](https://img.shields.io/github/downloads/heitzmann/gdspy/total.svg)](https://github.com/heitzmann/gdspy/releases)

Gdspy is a Python module for creation and manipulation of GDSII stream files.
Key features for the creation of complex CAD layouts are included:

* Boolean operations on polygons (AND, OR, NOT, XOR) based on clipping algorithm
* Polygon offset (inward and outward rescaling of polygons)
* Efficient point-in-polygon solutions for large array sets

Gdspy also includes a simple layout viewer.

Typical applications of Gdspy are in the fields of electronic chip design, planar lightwave circuit design, and mechanical engineering.


## Future of Gdspy

In trying to improve the performance of Gdspy for large layouts, we ended up concluding that the best way to reach our goal was to rewrite the critical parts of the library as a C extension.
It turns out that beside obvious functions, method calling has a big impact in performance due to the overhead it introduces.
The best solution was to re-design the whole project as a C++ library with a thin Python wrapper: thus was born [Gdstk, the GDSII Tool Kit](https://github.com/heitzmann/gdstk).

Therefore, version 1.6 will be the last major release of Gdspy, with development focused only on bug fixes.
Users are encouraged to move from Gdspy to Gdstk: although their API is not 100% compatible, the new module should be familiar enough to allow a quick transition.


## Installation

### Dependencies:

* [Python](https://www.python.org/) (tested with versions 2.7, 3.6, 3.7, and 3.8)
* [Numpy](http://numpy.scipy.org/)
* C compiler (needed only if built from source)
* Tkinter (optional: needed for the `LayoutViewer` GUI)
* [Sphinx](https://www.sphinx-doc.org/) (optional: to build the documentation)

### Linux / OS X

Option 1: using [pip](https://docs.python.org/3/installing/):

```sh
python -m pip install --user gdspy
```

Option 2: download the source from [github](https://github.com/heitzmann/gdspy) and build/install with:

```sh
python setup.py install
```

### Windows

The preferred option is to install pre-compiled binaries from [here](https://github.com/heitzmann/gdspy/releases).

Installation via `pip` and building from source as above are also possible, but an appropriate [build environment](https://wiki.python.org/moin/WindowsCompilers) is required for compilation of the C extension modules.


## Documentation

The complete documentation is available [here](http://gdspy.readthedocs.io/).

The source files can be found in the `docs` directory.


## Support

Help support Gdspy development by [donating via PayPal](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=JD2EUE2WPPBQQ)


## History of changes

### Version 1.6.12 (Jun 2, 2022)
* Fix in `Cell.get_texttypes`.
* Allow labels to inherit transforms through `get_labels`.

### Version 1.6.11 (Jan 14, 2022)
* Fix in `Cell.write_svg` when missing references.
* Speed improvements in `Cell.remove_polygons` (thanks Troy for the contribution).

### Version 1.6.10 (Nov 14, 2021)
* Fix in `Cell.get_polygons`

### Version 1.6.9 (Sep 23, 2021)
* Fix in `Cell.get_polygons` with specified layer and datatype.
* Raise error for duplicate cells when reading a GDSII file.

### Version 1.6.8 (Aug 2, 2021)
* Fix in `boolean` for complex geometries that freeze the operation.

### Version 1.6.7 (Jul 14, 2021)
* Fixes in `boolean` for bugs with self-intersecting holes and holes horizontal edges.
* Fix bug in warning message.

### Version 1.6.6 (Jun 09, 2021)
* Fix error in `Path.smooth` not finding `_hobby` function.
* Allow precision specification in SVG output.

### Version 1.6.5 (Jun 08, 2021)
* Support GDSII files with 0-padding at the end.
* Allow fixing and modifying GDSII file timestamps.
* *Thanks Troy Tamas and Joaquin Matres for the fixes*

### Version 1.6.4 (Apr 23, 2021)
* Fix missing module import (thanks Troy Tamas for the fix).

### Version 1.6.3 (Dec 28, 2020)
* Fix bounding box edge case (thanks Troy Tamas for the fix).

### Version 1.6.2 (Dec 18, 2020)
* More efficient bounding box calculation (thanks to Troy Tamas for the contribution).
* Fix Label creation bug.

### Version 1.6.1 (Oct 22, 2020)
* Fix SVG output when `Label` contains special characters.

### Version 1.6 (Aug 12, 2020)
* Added support for element properties.
* Added transformation support to `Cell.copy`.
* Layer/datatype filtering in `get_polygons` for `Cell`, `CellReference` and `CellArray`.
* Layer/datatype filtering in `LayoutViewer`.
* Removed global cache `_bounding_boxes`.  Only cells cache their bounding boxes.
* Bug fixes (thanks Daniel Hwang for the contributions).
* Bug fix in `Cell.copy` where the whole dependency tree would be copied on a deep copy creation.

### Version 1.5.2 (Feb 01, 2020)
* Added support for importing GDSII files containing BOX elements.
* Bug fix in `GdsLibrary.extract` (thanks collineps for finding the problem).

### Version 1.5 (Dec 20, 2019)
* New `Cell.write_svg` function to export an SVG image of the cell.
* New `GdsLibrary.new_cell` function to quickly create and add cells to a library.
* `GdsLibrary.add` can update references when a cell is overwritten.
* Added `GdsLibrary.remove` to allow cells to be properly removed from libraries.
* Added `GdsLibrary.rename_cell` to rename cells in libraries.
* Added `GdsLibrary.replace_references` to easily replace referenced cells in libraries.
* `GdsLibrary.add` can add dependencies recursively.
* Iterating over `GdsLibrary` objects yields all its cells.
* Iterating over `Cell` objects yield all its polygons, paths, labels and references.
* Breaking change to `*.to_gds` functions in order to improve write efficiency (this should not be a problem for most users, since `gdspy.write_gds` and `Cell.write_gds` remain the same).
* Breaking change: renamed `GdsLibrary.cell_dict` to `GdsLibrary.cells`.
* Deprecated: `gdspy.current_library`, `gdspy.write_gds`, `gdspy.fast_boolen`, `GdsLibrary.extract`.
* Bug fixes and better tests for `FlexPath` and `RobustPath`.

### Version 1.4.3 (Nov 11, 2019)
* Bug fix for `FlexPath` and `RobustPath` references.

### Version 1.4.2 (Oct 01, 2019)
* Bug fix in `FlexPath`.

### Version 1.4.1 (Sep 20, 2019)
* Bug fixes (thanks to DerekK88 and Sequencer for the patches).

### Version 1.4 (May 18, 2019)
* Revised [documentation](http://gdspy.readthedocs.io/).
* New `FlexPath` and `RobustPath` classes: more efficient path generation when using the original GDSII path specification.
* New `Curve` class: SVG-like polygon creation.
* Added `PolygonSet.mirror` (thanks to Daan Waardenburg for the contribution).
* Added `Path.bezier` to create paths based on Bézier curves.
* Added `Path.smooth` to create paths based on smooth interpolating curves.
* Added `get_gds_units` to get units used in a GDSII file without loading.
* Added `get_binary_cells` to load only the binary GDSII representation of cell from a file.
* Added argument `tolerance` to `Round`, `Path.arc`, `Path.turn`, and `Path.parametric` to automatically control the number of points in the final polygons.
* Added argument `binary_cells` to GDSII writing functions to support `get_binary_cells`.
* Added argument `rename_template` to `GdsLibrary.read_gds` for flexible cell renaming (thanks to @yoshi74ls181 for the contribution).
* Changed return value of `slice` to avoid creating empty `PolygonSet`.
* Added argument `timestamp` to GDSII writing functions.
* Improved `Round` to support creating ellipses.
* Added support for unlimited number of points per polygon.
* Added support for BGNEXTN and ENDEXTN when reading a GDSII file.
* Polygon creation warnings are now controlled by `poly_warnings`.
* Incorrect `anchor` in `Label` now raises an error, instead of emitting a warning.
* Added correct support for radius in `PolygonSet.fillet` on a per-vertex basis.
* Speed improvements in GDSII file generation (thanks to @fbeutel for the contribution) and geometry creation.
* Font rendering example using [matplotlib](https://matplotlib.org/) (thanks Hernan Pastoriza for the contribution).
* Expanded test suite.

### Version 1.3.2 (Mar 14, 2019)
* Small fix for building on Mac OS X Mojave.

### Version 1.3.1 (Jun 29, 2018)
* `PolygonSet` becomes the base class for all polygons, in particular `Polygon` and `Rectangle`.
* Added `Cell.remove_polygons` and `Cell.remove_labels` functions to allow filtering a cell contents based, for example, on each element's layer.
* Added `PolygonSet.scale` utility method.
* Added `PolygonSet.get_bounding_box` utility method.
* Added argument `timestamp` to `Cell.to_gds`, `GdsLibrary.write_gds` and `GdsWriter`.
* Added `unit` and `precision` arguments to `GdsLibrary` initialization and removed from its `write_gds` method.
* Changed the meaning of argument `unit` in `GdsLibrary.read_gds`.
* Improved `slice` to avoid errors when slicing in multiple positions at once.
* Improved `PolygonSet.fracture` to reduce number of function calls.
* Removed incorrect absolute flags for magnification and rotation in `CellReference` and `CellArray`.
* Minor bug fixes.
* Documentation fixes.
* Removed deprecated classes and functions.

### Version 1.2.1 (Dec 5, 2017)
* `GdsLibrary` can be created directly from a GDSII file
* Added return value to `GdsLibrary.read_gds`
* Fixed return value of `GdsLibrary.add`

### Version 1.2 (Oct 21, 2017)
* Added new `gdsii_hash` function.
* Added `precision` parameter to `_chop`, `Polygon.fracture`, `Polygon.fillet`, `PolygonSet.fracture`, `PolygonSet.fillet`, and `slice`.
* Included labels in flatten operations (added `get_labels` to `Cell`, `CellReference`, and `CellArray`).
* Fixed bug in the bounding box cache of reference copies.
* Fixed bug in `_chop` that affected `Polygon.fracture`, `PolygonSet.fracture`, and `slice`.
* Other minor bug fixes.

### Version 1.1.2 (Mar 19, 2017)
* Update clipper library to 6.4.2 to fix bugs introduced in the last update.
* License change to Boost Software License v1.0.

### Version 1.1.1 (Jan 27, 2017)
* Patch to fix installation issue (missing README file in zip).

### Version 1.1 (Jan 20, 2017)
* Introduction of `GdsLibrary` to allow user to work with multiple library simultaneously.
* Deprecated `GdsImport` in favor of `GdsLibrary`.
* Renamed `gds_print` to `write_gds` and `GdsPrint` to `GdsWriter`.
* Development changed to Python 3 (Python 2 supported via [python-future](http://python-future.org/)).
* Added photonics example.
* Added test suite.
* Clipper library updated to last version.
* Fixed `inside` function sometimes reversing the order of the output.
* Fixed rounding error in `fast_boolean`.
* Fixed argument `deep_copy` being inverted in `Cell.copy`.
* Bug fixes introduced by numpy (thanks to Adam McCaughan for the contribution).

### Version 1.0 (Sep 11, 2016)
* Changed to "new style" classes (thanks to Adam McCaughan for the contribution).
* Added	a per-point radius specification for `Polygon.fillet` (thanks to Adam McCaughan for the contribution).
* Added `inside` fucntion to perform point-in-polygon tests (thanks to @okianus for the contribution).
* Moved from distutils to setuptools for better Windows support.

### Version 0.9 (Jul 17, 2016)
* Added option to join polygons before applying an `offset`.
* Added a `translate` method to geometric entities (thanks John Bell for the commit).
* Bug fixes.

### Version 0.8.1 (May 6, 2016)
* New `fast_boolean` function based on the [Clipper](http://www.angusj.com/delphi/clipper.php) library with much better performance than the old `boolean`.
* Changed `offset` signature to also use the [Clipper](http://www.angusj.com/delphi/clipper.php) library (this change **breaks compatibility** with previous versions).
* Bug fix for error when importing some labels from GDSII files.

### Version 0.7.1 (June 26, 2015)
* Rebased to GitHub.
* Changed source structure and documentation.

### Version 0.7 (June 12, 2015)
* New feature: `offset` function.
* New `GdsPrint` class for incremental GDSII creation (thanks to Jack Sankey for the contribution).

### Version 0.6 (March 31, 2014)
* Default number of points for `Round`, `Path.arc`, and `Path.turn` changed to resolution of 0.01 drawing units.
* `Path.parametric` accepts callable `final_distance` and `final_width` for non-linear tapering.
* Added argument `ends` to `PolyPath`.
* Added (limited) support for PATHTYPE in `GdsImport`.
* A warning is issued when a `Path` curve has width larger than twice its radius (self-intersecting polygon).
* Added a random offset to the patterns in `LayoutViewer`.
* `LayoutViewer` shows cell labels for referenced cells.
* `get_polygons` returns (referenced) cell name if `depth` < 1 and `by_spec` is True.
* Bug fix in `get_bounding_box` when empty cells are referenced.
* Bug fixes in `GdsImport` and many speed improvements in bounding box calculations (thanks to Gene Hilton for the patch).

### Version 0.5 (October 30, 2013) - NOT COMPATIBLE WITH PREVIOUS VERSIONS
* Major `LayoutViewer` improvements (not backwards compatible).
* The layer argument has been repositioned in the argument list in all functions (not backwards compatible).  
* Renamed argument `by_layer` to `by_spec` (not backwards compatible).
* Error is raised for polygons with more vertices than possible in the GDSII format.
* Removed the global state variable for default datatype.
* Added `get_datatypes` to `Cell`.
* Added argument `single_datatype` to `Cell.flatten`.
* Removed `gds_image` and dropped the optional PIL dependency.

### Version 0.4.1 (June 5, 2013)
* Added argument `axis_offset` to `Path.segment` allowing creation of asymmetric tapers.
* Added missing argument  `x_reflection` to `Label`.
* Created a global state variable to override the default datatype.
* Bug fix in `CellArray.get_bounding_box` (thanks to George McLean for the fix)

### Version 0.4 (October 25, 2012)
* `Cell.get_bounding_box` returns `None` for empty cells.
* Added a cache for bounding boxes for faster computation, especially for references.
* Added support for text elements with `Label` class.
* Improved the emission of warnings.
* Added a tolerance parameter to `boolean`.
* Added better print descriptions to classes.
* Bug fixes in boolean involving results with multiple holes.

### Version 0.3.1 (May 24, 2012)
* Bug fix in the fracture method for `PolygonSet`.

### Version 0.3a (May 03, 2012)
* Bug fix in the fracture method for `Polygon` and `PolygonSet`.

### Version 0.3 (April 25, 2012)
* Support for Python 3.2 and 2.7
* Further improvements to the `boolean` function via caching.
* Added methods `get_bounding_box` and `get_layers` to `Cell`.
* Added method `top_level` to `GdsImport`.
* Added support for importing GDSII path elements.
* Added an argument to control the verbosity of the import function.
* Layer -1 (referenced cells) sent to the bottom of the layer list by default in `LayoutViewer`
* The text and background of the layer list in `LayoutViewer` now reflect the colors of the outlines and canvas backgroung.
* Changed default background color in `LayoutViewer`
* Thanks to Gene Hilton for the contributions!

### Version 0.2.9 (December 14, 2011)
* Attribute `Cell.cell_list` changed to `Cell.cell_dict`.
* Changed the signature of the operation in `boolean`.
* Order of cells passed to `LayoutViewer` is now respected in the GUI.
* Complete re-implementation of the boolean function as a C extension for improved performance.
* Removed precision argument in `boolean`. It is fixed at 1e-13 for merging close points, otherwise machine precision is used.
* `gds_image` now accepts cell names as input.
* Added optional argument `depth` to `get_polygons`
* Added option to convert layers and datatypes in imported GDSII cells.
* Argument `exclude_layers` from `LayoutViewer` changed to `hidden_layers` and behavior changed accordingly. 
* Shift + Right-clicking on a layer the layer-list of `LayoutVIewer` hides/unhides all other layers.
* New buttons to zoom in and out in `LayoutViewer`.
* Referenced cells below a configurable depth are now represented by theirs bounding boxes in `LayoutViewer`.

### Version 0.2.8 (June 21, 2011)
* GDSII file import
* GDSII output automatically include required referenced cells.
* `gds_print` also accepts file name as input.
* Outlines are visible by default in `LayoutViewer`.
* Added background color option in `LayoutViewer`.
* Right-clicking on the layer list hides/unhides the target layer in `LayoutViewer`.
* `Cell.cell_list` is now a dictionary indexed by name, instead of a list.
* Added option to exclude created cells from the global list of cells kept in `Cell.cell_list`.
* `CellReference` and `CellArray` accept name of cells as input.
* Submodules lost their own `__version__`.

### Version 0.2.7 (April 2, 2011)
* Bug fixed in the `boolean`, which affected the way polygons with more vertices then the maximum were fractured.
* `gds_image` accepts an extra color argument for the image background.
* Screenshots takes from `LayoutViewer` have the same background color as the viewer.
* The functions `boolean` and `slice` now also accept `CellReference` and `CellArray` as input.
* Added the method `fracture` to `Polygon` and `PolygonSet` to automatically slice polygons into parts with a predefined maximal number of vertices.
* Added the method `fillet` to `Polygon` and `PolygonSet` to round corners of polygons.

### Version 0.2.6 (February 28, 2011)
* When saving a GDSII file, `ValueError` is raised if cell names are duplicated.
* Save screenshot from `LayoutViewer`.
* `gds_image` accepts cells, instead of lists.
* Outlines supported by `gds_image`.
* `LayoutViewer` stores bounding box information for all visited layers to save rendering time.

### Version 0.2.5 (December 10, 2010)
* Empty cells no longer break the LayoutViewer.
* Removed the `gds_view` function, superseded by the LayoutViewer, along with all dependencies to matplotlib.
* Fixed a bug in `boolean` which affected polygons with series of collinear vertices.
* Added a function to `slice` polygons along straight lines parallel to an axis.

### Version 0.2.4 (September 04, 2010)
* Added shortcut to Extents in LayoutViewer: `Home` or `a` keys.
* `PolygonSet` is the new base class for `Round`, which might bring some incompatibility issues with older scripts.
* `Round` elements, `PolyPath`, `L1Path`, and `Path arc`, `turn` and `parametric` sections are now automatically fractured into pieces defined by a maximal number of points.
* Default value for `max_points` in boolean changed to 199.
* Removed the flag to disable the warning about polygons with more than 199 vertices.  The warning is shown only for `Polygon` and `PolygonSet`.
* Fixed a bug impeding parallel `parametric` paths to change their distance to each other.

### Version 0.2.3 (August 09, 2010)
* Added the `PolyPath` class to easily create paths with sharp corners.
* Allow `None` as item in the colors parameter of `LayoutViewer` to make layers invisible.
* Added color outline mode to `LayoutViewer` (change outline color with the shift key pressed)
* Increased the scroll region of the `LayoutViewer` canvas
* Added a fast scroll mode: control + drag 2nd mouse button
* Created a new sample script

### Version 0.2.2 (July 29, 2010)
* Changed the cursor inside `LayoutViewer` to standard arrow.
* Fixed bugs with the windows version of `LayoutViewer` (mouse wheel and ruler tool).

### Version 0.2.1 (July 29, 2010)
* Bug fix: `gds_image` displays an error message instead of crashing when `PIL` is not found.
* Added class `LayoutViewer`, which uses Tkinter (included in all Python distributions) to display the GDSII layout with better controls then the `gds_view` function. This eliminates the `matplotlib` requirement for the viewer functionality.
* New layer colors extending layers 0 to 63.

### Version 0.2.0 (July 19, 2010)
* Fixed a bug on the `turn` method of `Path`.
* Fixed a bug on the `boolean` function that would give an error when not using `Polygon` or `PolygonSet` as input objects.
* Added the method `get_polygons` to `Cell`, `CellReference` and `CellArray`.
* Added a copy method to `Cell`.
* Added a `flatten` method to `Cell` to remove references (or array references) to other cells.
* Fracture `boolean` output polygons based on the number of vertices to respect the 199 GDSII limit.  

### Version 0.1.9 (June 04, 2010)
* Added `L1Path` class for Manhattan geometry (L1 norm) paths.

### Version 0.1.8 (May 10, 2010)
* Removed the argument `fill` from `gds_view` and added a more flexible one: `style`.
* Fixed a rounding error on the `boolean` operator affecting polygons with holes.
* Added a rotate method to `PolygonSet`.
* Added a warning when `PolygonSet` has more than 199 points
* Added a flag to disable the warning about polygons with more than 199 points.
* Added a `turn` method to `Path`, which is easier to use than `arc`.
* Added a direction attribute to `Path` to keep the information used by the `segment` and `turn` methods.

### Version 0.1.7 (April 12, 2010)
* New visualization option: save the geometry directly to an image file (lower memory use).
* New functionality added: boolean operations on polygons (polygon clipping).
* All classes were adapted to work with the boolean operations.
* The attribute size in the initializer of class `Text` does not have a default value any longer.
* The name of the argument `format` in the function `gds_view` was changed to `fill` (to avoid confusion with the built-in function `format`).

### Version 0.1.6 (December 15,  2009)
* Sample script now include comments and creates an easier to understand GDSII example.
* Improved floating point to integer rounding, which fixes the unit errors at the last digit of the precision in the GDSII file.
* Fixed the font for character 5.
* Added a flag to `gds_view` to avoid the automatic call to `matplotlib.pyplot.show()`.
* In `gds_view`, if a layer number is greater than the number of formats defined, the formats are cycled.

### Version 0.1.5a (November 15, 2009)
* Class Text correctly interprets `\n` and `\t` characters.
* Better documentation format, using the Sphinx engine and the numpy format.

### Version 0.1.4 (October 5, 2009)
* Class `Text` re-written with a different font with no overlaps and correct size.

### Version 0.1.3a (July 29 2009)
* Fixed the function `to_gds` of class `Rectangle`.

### Version 0.1.3 (July 27, 2009)
* Added the datatype field to all elements of the GDSII structure.

### Version 0.1.2 (July 11, 2009)
* Added the `gds_view` function to display the GDSII structure using the matplotlib module.
* Fixed a rotation bug in the CellArray class.
* Module published under the GNU General Public License (GPL)

### Version 0.1.1 (May 12, 2009)
* Added attribute `cell_list` to class Cell to hold a list of all Cell created.
* Set the default argument `cells=Cell.cell_list` in the function `gds_print`.
* Added member to calculate the area for each element type.
* Added member to calculate the total area of a Cell or the area by layer.
* Included the possibility of creating objects in user-defined units, not only nanometers.

### Version 0.1.0 (May 1, 2009)
* Initial release.
