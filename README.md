# GDSPY README

[![Boost Software License - Version 1.0](https://img.shields.io/github/license/heitzmann/gdspy.svg)](http://www.boost.org/LICENSE_1_0.txt)
[![Documentation Status](https://readthedocs.org/projects/gdspy/badge/?version=latest)](http://gdspy.readthedocs.io/en/latest/?badge=latest)

Gdspy is a Python  module for creating/importing/merging GDSII stream files.
It includes key libraries for creating complex CAD layouts:

* Boolean operations on polygons (AND, OR, NOT, XOR) based on clipping algorithm
* Polygon offset (inward and outward rescaling of polygons)
* Efficient point-in-polygon solutions for large array sets

Gdspy also includes a simple layout viewer.

Typical applications of gdspy are in the fields of electronic chip design, planar lightwave circuit design, and mechanical engineering.

## Installation

### Dependencies:

* [Python](http://www.python.org/) (tested with versions 2.7, 3.5 and 3.6)
* [Numpy](http://numpy.scipy.org/)
* [Python-future](http://python-future.org/) (only for Python 2)
* C compiler (needed only if built from source)

### Linux / OS X

Option 1: using [pip](https://docs.python.org/3/installing/):

```sh
pip install gdspy
```

Option 2: download the source from [github](https://github.com/heitzmann/gdspy) and build/install with:

```sh
python setup.py install
```

### Windows

The preferred option is to install pre-compiled binaries from [here](https://github.com/heitzmann/gdspy/releases) for 32 and 64-bit systems.

Installation via `pip` and building from source as above are also possible, but an appropriate [build environment](https://wiki.python.org/moin/WindowsCompilers) is required for compilation of the C extension modules.

## Usage

The file [tutorial.py](https://github.com/heitzmann/gdspy/blob/master/examples/tutorial.py) in the `example` folder is a sample script to show the features provided by this module.

The complete module reference can be built from the sources in the `docs` folder with [Sphinx](http://sphinx-doc.org/).
It is also available [on-line](http://gdspy.readthedocs.io/en/latest/)

## Support

Help support gdspy development by [donating via PayPal](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=JD2EUE2WPPBQQ)

## History of changes

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
* [![Documentation Status](https://readthedocs.org/projects/gdspy/badge/?version=latest)](http://gdspy.readthedocs.io/en/latest/?badge=latest)Removed the argument `fill` from `gds_view` and added a more flexible one: `style`.
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
