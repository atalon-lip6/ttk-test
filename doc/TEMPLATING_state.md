TTK Developer Experience
========================

Author: Alexandre Talon
Date: December 202"

This document summarises the state of a work: templatising the triangulation
classes. The code is already templatise according to the type of the
triangulation: once a triangulation is passed to a ttk base module, the module
in general has code ready and optimised for the triangulation type of its
argument. This is made possible by generating at compilation time code for each
type of triangulation. At runtime a switch is done when calling a routine from
the module, so as to call the optimised vesion of that routine for the type of
triangulation of the current object.

This work aims at refining this by templating any triangulation by its
cardinality, which is so far constant thrhoughout the lifetime of an
object (unless it is reassigned). We now have 4 (0 to 3) versions of
the ExplicitTriangulation class for instance. We did not apply the
templatisation to the Triangle container class and to the
AbstractTriangulation base class. This would be between impossible and very
tedious in practice. Also, this would normally not improve the performance
since in all calls to ttk base modules we should explicitely pass an object
of a final triangulation.


## Already implemented
 * Each class derived from AbstractTriangulation is templated

* The calls from the vtk layer to the base modules now use the templated
 versions of the triangulation to pass their triangulation arguments:
    * ttkVtkTemplateMacro, ttkTypeMacroT, and other macros have been modified
    to use the new triangulation objects.

* The ExplicitTriangulation class now overrides more methods from
 AbstractTriangulation: the ones which behave differently according to the
 dimensionality of the object, apart from the preconditioning ones (not
 significant for performance improvements).

## Still TODO
 * Adapt a few base functions if they declare, use or take in argument some
 specialised triangulation (i.e. not of type AbstractTriangulation).
    * This was begun (but not tested) for BarycentricSubdivision
    * Also adapt the vtk modules calling them (easy).

 * Override some methods from AbstractTriangulation in all the derived
 triangulation classes except ExplicitTriangulation (already done).

 * Run the tests.

 * Do some benchmarking.

 * Check rather during compilation time, for type of triangulation in things using progressive/approximate topoloty

 * Maybe template multires topology also (and approximate, progressive)
    * Approxtopo and progressivetopo if templatées => simpler? 

 * For modules using setupTriangulation from multires, maybe do the templating inside multires setupTriangulation (and a new function?) instead of before calling setupTriangulation (less code)

 * Triangulation : isExplicit, isImplicit, etc. ?
    => DiscreteMorseSandwich.h


## Performance results:
 * 4-5 order of magnitude faster if calling just one function multiple times:
 the templatisation and use of if constexpr () works.
 * Apparently ~10% less time with DiscreteGradient: building th gradient and
 computing glyphs. On viscousFingering.vtu, using 1 thread.
    * Measured on the user time, averaging twenty runs.
    * Also measured with Valgrind: ~13% less time for RequestData, 15% for the BuildGradient call, 22% for FillGradientGlyphs, and 39% for the calls to GetTriangleInCenter (the optimised function).

 
